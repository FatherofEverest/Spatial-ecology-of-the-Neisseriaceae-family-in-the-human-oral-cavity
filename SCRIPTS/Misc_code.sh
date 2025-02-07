# variables
projectID=P_0003_Neisseriaceae
mainDIR=/workspace/jmarkwelchlab/$projectID
DIR_ContigsDB=04_CONTIGS_DB
DIR_SinglePROF=06_SINGLE_PROFILE
DIR_MergedPROF=07_MERGED_PROFILE
CONTIGS_DB=$DIR_ContigsDB/${projectID}-contigs.db
SAMPLE_NAME=${projectID}_TD
MERGED_PROFILE=$DIR_MergedPROF/$SAMPLE_NAME
SINGLE_PROFILES_DB=$DIR_SinglePROF/TD_*/PROFILE.db


# merged profiling per sample
anvi-merge -c $CONTIGS_DB --enforce-hierarchical-clustering -o $MERGED_PROFILE -S $SAMPLE_NAME $SINGLE_PROFILES_DB


```{bash, eval=FALSE}

output_file="TD_single_profile_details.txt"

for sample in 06_SINGLE_PROFILE/TD_HC_HMP_*/PROFILE.db
do
    if [ -f "$sample" ]; then
            echo -e "$sample" >> "$output_file"
            # Redirect both stdout and stderr to the file
            anvi-show-misc-data -p "$sample" >> "$output_file" 2>&1
    fi
done


echo -e "sample\tgroup" > TD_single_profile_groups.txt 
cat $output_file | grep 'TD_HC_HMP_\|* DATA GROUP' | sed -n -e '/06_SINGLE_PROFILE\/TD_HC_HMP_/{s|.*06_SINGLE_PROFILE/\(TD_HC_HMP_[^/]*\)/.*|\1|;h}' \
       -e '/\* DATA GROUP/{s|.*DATA GROUP "\(.*\)" WITH.*|\1|;H;x;s/\n/\t/;p}' >> TD_single_profile_groups.txt 

cat TD_single_profile_groups.txt | wc -l


# get list of TD samples
cat samples_id-QC_IDs.txt | grep "TD_" | sort > TD_samples_id-QC_IDs.txt

# compare the first column of TD_single_profile_groups.txt with the list of TD samples in TD_samples_id-QC_IDs.txt to find which samples are missing from TD_single_profile_groups.txt

# extract the First Column from TD_single_profile_groups.txt:
cut -f1 TD_single_profile_groups.txt > extracted_TD_single_profile_groups.txt

# sort Both Lists:
sort extracted_TD_single_profile_groups.txt > sorted_TD_single_profile_groups.txt
sort TD_samples_id-QC_IDs.txt > sorted_TD_samples_id-QC_IDs.txt

# compare the Lists with comm:
comm -23 sorted_TD_samples_id-QC_IDs.txt sorted_TD_single_profile_groups.txt > missing_samples.txt

```


Looks like the single profiles just randomly stopped for not particular reason, as all five samples the run log ends abbruptly at the line that reports the Number of nucleotides. For example,

[18 Nov 23 19:46:11] Number of nucleotides .............................: 486,609,095   

Unfortunately it appears that we have to find 8 other samples that failed profiling. The following code repeats the above analysis for the remaining oral sites. This will provide a separate file for each oral site that lists the sample IDs that are missing the "default" group in their profile databases. We can then append this list to the liost we have for the 5 TD samples and the 1 BM sample and re-run the single profiles for those samples. 

```{bash, eval=FALSE}

for site in BM PP PT TH KG PB SA
do

  output_file="${site}_single_profile_details.txt"

  for sample in 06_SINGLE_PROFILE/${site}_HC_HMP_*/PROFILE.db
  do
      if [ -f "$sample" ]; then
            echo -e "$sample" >> "$output_file"
            # Redirect both stdout and stderr to the file
            anvi-show-misc-data -p "$sample" >> "$output_file" 2>&1
      fi
  done


  echo -e "sample\tgroup" > ${site}_single_profile_groups.txt 
  cat $output_file | grep "${site}_HC_HMP_\|* DATA GROUP" | sed -n -e "/06_SINGLE_PROFILE\/${site}_HC_HMP_/{s|.*06_SINGLE_PROFILE/\(${site}_HC_HMP_[^/]*\)/.*|\1|;h}" \
       -e '/\* DATA GROUP/{s|.*DATA GROUP "\(.*\)" WITH.*|\1|;H;x;s/\n/\t/;p}' >> ${site}_single_profile_groups.txt 

  # get list of TD samples
  cat samples_id-QC_IDs.txt | grep "${site}_" | sort > ${site}_samples_id-QC_IDs.txt

  # extract the First Column from TD_single_profile_groups.txt:
  cut -f1 ${site}_single_profile_groups.txt > extracted_${site}_single_profile_groups.txt

  # sort Both Lists:
  sort extracted_${site}_single_profile_groups.txt > sorted_${site}_single_profile_groups.txt
  sort ${site}_samples_id-QC_IDs.txt > sorted_${site}_samples_id-QC_IDs.txt

  # compare the Lists with comm:
  comm -23 sorted_${site}_samples_id-QC_IDs.txt sorted_${site}_single_profile_groups.txt > ${site}_missing_samples.txt

done

```

