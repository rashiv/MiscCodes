#!/usr/bin/env bash

ls -la |grep ^d| awk '{if(NR>2) print $9}' > novels
# because first 2 lines have current and one up directories
for i in $(<novels)
do
  cd $i
  ls -la |grep ^d| awk '{if(NR>2) print $9}'> disks
  for j in $(<disks)
  do
    cd $j
    for file in *.wma
    do 
      ffmpeg -i "${file}"  -acodec libmp3lame -ab 192k "${file/.wma/.mp3}"
    done
    cd ..
    mp3wrap ${j}.mp3 $j/*.mp3
    mv ${j}_MP3WRAP.mp3 ${j}.mp3
# concatenating individual mp3 files in the disk
  done
  rm -f disks
  cd ..
done 
rm -f novels
