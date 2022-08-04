#!/bin/bash

while getopts 'i:' options
do
  case $options in
    i) fasta=$OPTARG ;;
  esac
done

# index the fasta input and get the sequence length
samtools faidx ${fasta}
seq_len=$(cut -f2 ${fasta}.fai)

mkdir -p _tmp_standardize 

# get the header 
head -1 ${fasta} > _tmp_standardize/header

## blast the fasta to itself to get the repeat alignments
## ideally, the IR sequences in plastome should be identical, and hence a perc_identity of 99 is used
## self match is removed and sorted by alignment length
blastn -query ${fasta} -subject ${fasta} -perc_identity 99 -evalue 0.00001 -outfmt '6 qseqid qstart qend sseqid sstart send length pident'|awk '$2!=$5&&$3!=$6 {print}'|sort  -k7,7nr -k2,2n  > _tmp_standardize/_tmp.blast.out
  
head -2 _tmp_standardize/_tmp.blast.out > _tmp_standardize/_tmp.blast.out2 
start=$(head -1 _tmp_standardize/_tmp.blast.out|awk '{print $2}')
end=$(head -1 _tmp_standardize/_tmp.blast.out|awk '{print $5}')
if [ $start -eq 1 ];then  
  awk '$3=='$seq_len' && $6-1=='$end' {print}' _tmp_standardize/_tmp.blast.out >> _tmp_standardize/_tmp.blast.out2 
  awk '$5=='$seq_len' && $2-1=='$end' {print}' _tmp_standardize/_tmp.blast.out >> _tmp_standardize/_tmp.blast.out2 
elif [ $end -eq $seq_len ];then
  awk '$2==1 && $5+1=='$start' {print}' _tmp_standardize/_tmp.blast.out >> _tmp_standardize/_tmp.blast.out2
  awk '$6==1 && $3+1=='$start' {print}' _tmp_standardize/_tmp.blast.out >> _tmp_standardize/_tmp.blast.out2
fi
mv _tmp_standardize/_tmp.blast.out2 _tmp_standardize/_tmp.blast.out

# if the fasta is split at inverted repeats, it would give more than 2 alignments (4 to be precise)
if [ $(grep -c "^" _tmp_standardize/_tmp.blast.out) -gt 2 ];then 

    echo "$(date) : Standardizing the input fasta file"
    ir_tmp=$(cat _tmp_standardize/_tmp.blast.out|awk '$3=='$seq_len' {print $6-1}')
    sc1=$(cat _tmp_standardize/_tmp.blast.out|awk '$5=='$ir_tmp' {print $6-$3}')
    sc2=$(cat _tmp_standardize/_tmp.blast.out|awk '$3=='$seq_len' {print $2-$5}')

    if [ $sc2 -gt $sc1 ];then 
        # do this if the fasta is split at IRa region 
        lsc=$(cat _tmp_standardize/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$5+1"-"$2-1}')
        ira_start=$(cat _tmp_standardize/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$2"-"$3}')
        ira_end=$(cat _tmp_standardize/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$2"-"$3}') 
        ssc=$(cat _tmp_standardize/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$3+1"-"$6-1}')
        irb_start=$(cat _tmp_standardize/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$6"-"$5}') 
        irb_end=$(cat _tmp_standardize/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$6"-"$5}')

    else 
        # do this if the fasta is split at IRb region 
        lsc=$(cat _tmp_standardize/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$3+1"-"$6-1}')
        ira_start=$(cat _tmp_standardize/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$6"-"$5}') 
        ira_end=$(cat _tmp_standardize/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$6"-"$5}')
        ssc=$(cat _tmp_standardize/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$5+1"-"$2-1}')
        irb_start=$(cat _tmp_standardize/_tmp.blast.out|awk '$3=='$seq_len' {print $1":"$2"-"$3}')
        irb_end=$(cat _tmp_standardize/_tmp.blast.out|awk '$5=='$ir_tmp' {print $1":"$2"-"$3}') 
    fi 

    samtools faidx ${fasta} "$irb_start" |sed "s/$irb_start/IRb/" > _tmp_standardize/irbs 
    samtools faidx ${fasta} "$irb_end" |sed 1d > _tmp_standardize/irbe  
    samtools faidx ${fasta} "$ira_start" |sed "s/$ira_start/IRa/" > _tmp_standardize/iras 
    samtools faidx ${fasta} "$ira_end" |sed 1d > _tmp_standardize/irae  

    cat _tmp_standardize/irbs _tmp_standardize/irbe > _tmp_standardize/IRb.fa 
    cat _tmp_standardize/iras _tmp_standardize/irae > _tmp_standardize/IRa.fa 

      # double checking the IR alignments 
      if [ "$(blastn -query _tmp_standardize/IRa.fa -subject _tmp_standardize/IRb.fa -perc_identity 99 -evalue 0.00001 -outfmt '6 qseqid qstart qend sseqid sstart send length pident' |awk 'NR==1&&$2==$6&&$3==$5 {print "Success"}')" == "Success" ] ;then 
        samtools faidx ${fasta} "$lsc" |sed 1d > _tmp_standardize/LSC.fa 
        samtools faidx ${fasta} "$ssc" |sed 1d > _tmp_standardize/SSC.fa 
        cat _tmp_standardize/IRa.fa | sed 1d > _tmp_standardize/_tmp.IRa.fa 
        cat _tmp_standardize/IRb.fa | sed 1d > _tmp_standardize/_tmp.IRb.fa
        cat _tmp_standardize/header _tmp_standardize/LSC.fa _tmp_standardize/_tmp.IRa.fa _tmp_standardize/SSC.fa _tmp_standardize/_tmp.IRb.fa | sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' > $(echo "${fasta%.*}").standardized.fa 
        rm -r _tmp_standardize/ ${fasta}.fai 
        echo "$(date) : Finished"
      fi

else 
    # if the end or start the sequence is there in the alignemnts, that means the IR locations needs to be adjusted
    if [ "$(cat _tmp_standardize/_tmp.blast.out| awk '$3=='$seq_len'||$2==1 {print "yes"}')" == "yes" ];then 
      if [ "$(cat _tmp_standardize/_tmp.blast.out| awk '$3=='$seq_len' {print "yes"}')" == "yes" ];then
        if [ $(awk '$3=='$seq_len' {print $6}' _tmp_standardize/_tmp.blast.out) -gt 50000 ];then
          echo "$(date) : Standardizing the input fasta file"
          sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D' ${fasta} > $(echo "${fasta%.*}").standardized.fa 
          rm -r _tmp_standardize/ ${fasta}.fai
          echo "$(date) : Finished"
        else 
          echo "$(date) : Standardizing the input fasta file"
          ira=$(cat _tmp_standardize/_tmp.blast.out |awk '$3=='$seq_len' {print $1":"$2"-"$3}')
          irb=$(cat _tmp_standardize/_tmp.blast.out |awk '$3=='$seq_len' {print $1":"$6"-"$5}')
          lsc=$(cat _tmp_standardize/_tmp.blast.out |awk '$3=='$seq_len' {print $1":"$5+1"-"$2-1}')  
          ssc=$(cat _tmp_standardize/_tmp.blast.out |awk '$3=='$seq_len' {print $1":1-"$6-1}')  

          samtools faidx ${fasta} "$irb" |sed 1d > _tmp_standardize/IRb.fa
          samtools faidx ${fasta} "$ira" |sed 1d > _tmp_standardize/IRa.fa 
          samtools faidx ${fasta} "$ssc" |sed 1d > _tmp_standardize/ssc.fa 
          samtools faidx ${fasta} "$lsc" |sed 1d > _tmp_standardize/lsc.fa 

          cat _tmp_standardize/header _tmp_standardize/lsc.fa _tmp_standardize/IRa.fa _tmp_standardize/ssc.fa _tmp_standardize/IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'> $(echo "${fasta%.*}").standardized.fa 
          rm -r _tmp_standardize/ ${fasta}.fai
          echo "$(date) : Finished"
        fi 
      fi 
      if [ "$(cat _tmp_standardize/_tmp.blast.out| awk '$2==1 {print "yes"}')" == "yes" ];then
        if [ $(cat _tmp_standardize/_tmp.blast.out |awk '$2==1 {print '$seq_len'-$5}') -gt 50000 ];then 
          echo "$(date) : Standardizing the input fasta file"
          ira=$(cat _tmp_standardize/_tmp.blast.out |awk '$2==1 {print $1":"$2"-"$3}')
          irb=$(cat _tmp_standardize/_tmp.blast.out |awk '$2==1 {print $1":"$6"-"$5}')
          lsc=$(cat _tmp_standardize/_tmp.blast.out |awk '$2==1 {print $1":"$5+1"-"'$seq_len'}')  
          ssc=$(cat _tmp_standardize/_tmp.blast.out |awk '$2==1 {print $1":"$3+1"-"$6-1}')  

          samtools faidx ${fasta} "$irb" |sed 1d > _tmp_standardize/IRb.fa
          samtools faidx ${fasta} "$ira" |sed 1d > _tmp_standardize/IRa.fa 
          samtools faidx ${fasta} "$ssc" |sed 1d > _tmp_standardize/ssc.fa 
          samtools faidx ${fasta} "$lsc" |sed 1d > _tmp_standardize/lsc.fa 

          cat _tmp_standardize/header _tmp_standardize/lsc.fa _tmp_standardize/IRa.fa _tmp_standardize/ssc.fa _tmp_standardize/IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'> $(echo "${fasta%.*}").standardized.fa 
          rm -r _tmp_standardize/ ${fasta}.fai
          echo "$(date) : Finished"
        else 
          echo "$(date) : Standardizing the input fasta file"
          irb=$(cat _tmp_standardize/_tmp.blast.out |awk '$2==1 {print $1":"$2"-"$3}')
          ira=$(cat _tmp_standardize/_tmp.blast.out |awk '$2==1 {print $1":"$6"-"$5}')
          ssc=$(cat _tmp_standardize/_tmp.blast.out |awk '$2==1 {print $1":"$5+1"-"'$seq_len'}')  
          lsc=$(cat _tmp_standardize/_tmp.blast.out |awk '$2==1 {print $1":"$3+1"-"$6-1}')  

          samtools faidx ${fasta} "$irb" |sed 1d > _tmp_standardize/IRb.fa
          samtools faidx ${fasta} "$ira" |sed 1d > _tmp_standardize/IRa.fa 
          samtools faidx ${fasta} "$ssc" |sed 1d > _tmp_standardize/ssc.fa 
          samtools faidx ${fasta} "$lsc" |sed 1d > _tmp_standardize/lsc.fa 

          cat _tmp_standardize/header _tmp_standardize/lsc.fa _tmp_standardize/IRa.fa _tmp_standardize/ssc.fa _tmp_standardize/IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'> $(echo "${fasta%.*}").standardized.fa 
          rm -r _tmp_standardize/ ${fasta}.fai
          echo "$(date) : Finished"
        fi 
      fi  
    else 
          # otherwise check if it is split at LSC or SSC
          echo "$(date) : Standardizing the input fasta file"
          end1=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print '$seq_len'-$3}')
          end2=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $6}')
          sum=$(($end1 + $end2))
          mid=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $2-$5}')

          if [ $mid -gt $sum ];then 
            # the fasta is split at SSC region
            lsc=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$5+1"-"$2-1}')
            ira=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$2"-"$3}')
            ssc_start=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$3+1"-"'$seq_len'}')
            ssc_end=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":1-"$6-1}')
            irb=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$6"-"$5}')

            samtools faidx ${fasta} "$irb" |sed 1d > _tmp_standardize/IRb.fa
            samtools faidx ${fasta} "$ira" |sed 1d > _tmp_standardize/IRa.fa 
            samtools faidx ${fasta} "$ssc_start" |sed 1d > _tmp_standardize/sscs.fa 
            samtools faidx ${fasta} "$ssc_end" |sed 1d > _tmp_standardize/ssce.fa 
            samtools faidx ${fasta} "$lsc" |sed 1d > _tmp_standardize/lsc.fa 

            cat _tmp_standardize/header _tmp_standardize/lsc.fa _tmp_standardize/IRa.fa _tmp_standardize/sscs.fa _tmp_standardize/ssce.fa _tmp_standardize/IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'> $(echo "${fasta%.*}").standardized.fa 
            rm -r _tmp_standardize/ ${fasta}.fai
            echo "$(date) : Finished"

          else 
                # the fasta is split at LSC region
                lsc_start=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$3+1"-"'$seq_len'}')
                lsc_end=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":1-"$6-1}')
                ira=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$6"-"$5}')
                ssc=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$5+1"-"$2-1}')
                irb=$(cat _tmp_standardize/_tmp.blast.out | awk '$2>$5&&$3>$6 {print $1":"$2"-"$3}')

                samtools faidx ${fasta} "$irb" |sed 1d > _tmp_standardize/IRb.fa
                samtools faidx ${fasta} "$ira" |sed 1d > _tmp_standardize/IRa.fa 
                samtools faidx ${fasta} "$lsc_start" |sed 1d > _tmp_standardize/lscs.fa 
                samtools faidx ${fasta} "$lsc_end" |sed 1d > _tmp_standardize/lsce.fa 
                samtools faidx ${fasta} "$ssc" |sed 1d > _tmp_standardize/ssc.fa 

                cat _tmp_standardize/header _tmp_standardize/lscs.fa _tmp_standardize/lsce.fa _tmp_standardize/IRa.fa _tmp_standardize/ssc.fa _tmp_standardize/IRb.fa |sed ':a; $!N; /^>/!s/\n\([^>]\)/\1/; ta; P; D'> $(echo "${fasta%.*}").standardized.fa 
                rm -r _tmp_standardize/ ${fasta}.fai
                echo "$(date) : Finished"
          fi
    fi
fi
