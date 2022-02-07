source '/cromwell_root/gcs_transfer.sh'

timestamped_message 'Delocalization script execution started...'

# encode-processing
delocalize_2b58939d5690e9e4dbcb374ad0b13fa0=(
  "encode-processing"       # project
  "3"   # max attempts
  "0" # parallel composite upload threshold, will not be used for directory types
  "file"
  "gs://encode-processing/caper_out_v04_05/segway/0105c048-0769-4499-9b37-d72d693a7470/call-interpretation/memory_retry_rc"
  "/cromwell_root/memory_retry_rc"
  "optional"
  "text/plain; charset=UTF-8"
  "file"
  "gs://encode-processing/caper_out_v04_05/segway/0105c048-0769-4499-9b37-d72d693a7470/call-interpretation/rc"
  "/cromwell_root/rc"
  "required"
  "text/plain; charset=UTF-8"
  "file"
  "gs://encode-processing/caper_out_v04_05/segway/0105c048-0769-4499-9b37-d72d693a7470/call-interpretation/monitoring.log"
  "/cromwell_root/monitoring.log"
  "required"
  "text/plain; charset=UTF-8"
  "file"
  "gs://encode-processing/caper_out_v04_05/segway/0105c048-0769-4499-9b37-d72d693a7470/call-interpretation/stdout"
  "/cromwell_root/stdout"
  "required"
  "text/plain; charset=UTF-8"
  "file"
  "gs://encode-processing/caper_out_v04_05/segway/0105c048-0769-4499-9b37-d72d693a7470/call-interpretation/stderr"
  "/cromwell_root/stderr"
  "required"
  "text/plain; charset=UTF-8"
  "file"
  "gs://encode-processing/caper_out_v04_05/segway/0105c048-0769-4499-9b37-d72d693a7470/call-interpretation/interpretation-output/classification/sample/probs.txt"
  "/cromwell_root/interpretation-output/classification/sample/probs.txt"
  "required"
  ""
  "file"
  "gs://encode-processing/caper_out_v04_05/segway/0105c048-0769-4499-9b37-d72d693a7470/call-interpretation/interpretation-output/classification/sample/classifier_data.tab"
  "/cromwell_root/interpretation-output/classification/sample/classifier_data.tab"
  "required"
  ""
  "file"
  "gs://encode-processing/caper_out_v04_05/segway/0105c048-0769-4499-9b37-d72d693a7470/call-interpretation/interpretation-output/classification/sample/mnemonics.txt"
  "/cromwell_root/interpretation-output/classification/sample/mnemonics.txt"
  "required"
  ""
  "file"
  "gs://encode-processing/caper_out_v04_05/segway/0105c048-0769-4499-9b37-d72d693a7470/call-interpretation/glob-4d10984b9a0d95111a61cf43c925df9e.list"
  "/cromwell_root/glob-4d10984b9a0d95111a61cf43c925df9e.list"
  "required"
  ""
  "directory"
  "gs://encode-processing/caper_out_v04_05/segway/0105c048-0769-4499-9b37-d72d693a7470/call-interpretation/glob-4d10984b9a0d95111a61cf43c925df9e/"
  "/cromwell_root/glob-4d10984b9a0d95111a61cf43c925df9e"
  "required"
  ""
)

delocalize "${delocalize_2b58939d5690e9e4dbcb374ad0b13fa0[@]}"
      
timestamped_message 'Delocalization script execution complete.'