#! bin/bash

SOURCE="/scratch/mepicard/physio_data/emotionvideos"
OUTDIR="/scratch/mepicard/physio_data/emotionvideos/processed"

#python neuromod_process.py $SOURCE "sub-02" "ses-001" $OUTDIR True
#python neuromod_process.py $SOURCE "sub-02" "ses-002" $OUTDIR True
#python neuromod_process.py $SOURCE "sub-02" "ses-003" $OUTDIR True
#python neuromod_process.py $SOURCE "sub-02" "ses-004" $OUTDIR True
#python neuromod_process.py $SOURCE "sub-02" "ses-005" $OUTDIR True
#python neuromod_process.py $SOURCE "sub-02" "ses-006" $OUTDIR True
#python neuromod_process.py $SOURCE "sub-02" "ses-007" $OUTDIR True

#python neuromod_process.py $SOURCE "sub-03" "ses-002" $OUTDIR True
#python neuromod_process.py $SOURCE "sub-03" "ses-003" $OUTDIR True
python neuromod_process.py $SOURCE "sub-03" "ses-004" $OUTDIR True
python neuromod_process.py $SOURCE "sub-03" "ses-005" $OUTDIR True
python neuromod_process.py $SOURCE "sub-03" "ses-006" $OUTDIR True
python neuromod_process.py $SOURCE "sub-03" "ses-007" $OUTDIR True

python neuromod_process.py $SOURCE "sub-05" "ses-001" $OUTDIR True
python neuromod_process.py $SOURCE "sub-05" "ses-002" $OUTDIR True
python neuromod_process.py $SOURCE "sub-05" "ses-003" $OUTDIR True
python neuromod_process.py $SOURCE "sub-05" "ses-004" $OUTDIR True




