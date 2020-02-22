# render script

set n_frames [molinfo top get numframes]

exec rm -rf images
exec mkdir -p images

display resize 800 800

# Rewind trajectory and record one more time
for { set frame 0 } { $frame < $n_frames } { incr frame } {
    animate goto $frame
    set d [format %05d $frame]
    render snapshot images/$d.tga
}

exec ffmpeg -y -framerate 25 -i images/%05d.tga -c:v libx264 -strict -2 -preset slow -pix_fmt yuv420p -vb 20M -s 800x800 tip3p.avi
