#!/bin/bash
#
# renders movie

tag=$1
interlace=$2

if [ "$1" == "" ]; then
  echo "usage: ./xcreate_ffmpeg_movie.sh [tag,e.g. marsUHD] [interlace,e.g. 6]"
  exit
fi

echo "start: `date`"

echo
echo "converting to movie..."
echo
rm -f movie.mp4

## ffmpeg options:
# quality options
# compression
comp=h264

# H.264 constant rate factor
crf=24  # better quality: 20

#pixel format
pxfmt=yuvj420p

# framerate
if [ "$interlace" == "" ]; then
  framerate=20
else
  framerate=$(( 20 * interlace ))
fi

## slow framerate
#framerate=15    # default DT = 0.003 and NSTEP per movie snapshot 100 -> time per frame = 0.003 * 100 = 0.3 s
#framerate=10    # since time per frame here is DT * NSTEP per movie snapshot            = 0.02 * 40   = 0.8 s -> slow down a bit with lower framerate
framerate=8      # time per snapshot 0.025 * 200 = 5 s

#if [ $framerate -gt 80 ]; then framerate=80; fi
if [ $framerate -gt 100 ]; then framerate=100; fi
#if [ $framerate -gt 120 ]; then framerate=120; fi  # framerates of 120 can trigger Quicktimes slow motion features...

echo "options: framerate = $framerate"
echo

# quality in case crf option not available
if [[ $tag =~ "UHD" ]];then
bitrate="-b:v 30000k"
else
bitrate="-b:v 9000k"
fi

# for jpg example:
# > ffmpeg -framerate 10 -f image2 -pattern_type glob -i 'movie/frame.*.jpg' -c:v h264 -crf 24 -pix_fmt yuv420p movie.mp4

echo "running:"
echo "> ffmpeg -framerate $framerate -f image2 -pattern_type glob -i 'movie/frame.*.jpg' -c:v $comp -crf $crf -pix_fmt $pxfmt movie.mp4"
echo

ffmpeg -framerate $framerate -f image2 -pattern_type glob -i 'movie/frame.*.jpg' -c:v $comp -crf $crf -pix_fmt $pxfmt movie.mp4

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# to get movie info
#ffmpeg -i movie.mp4

# to shorten movie time
# see: https://trac.ffmpeg.org/wiki/How%20to%20speed%20up%20/%20slow%20down%20a%20video
#ffmpeg -i movie.mp4 -filter:v "setpts=0.5*PTS" output.mp4

# to smooth movie
#ffmpeg -i movie.mp4 -filter:v "minterpolate='mi_mode=mci:mc_mode=aobmc:vsbmc=1:fps=120'" output.mp4

echo
echo "written to file: movie.mp4"
echo
echo
