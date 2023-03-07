#! /bin/sh

for FILE in $1*
do
tmp=${FILE##*/}
ffprobe -hide_banner -show_entries "stream=filename,duration,nb_frames,r_frame_rate" -of json $FILE > ${2}${tmp%.*}.json
done

