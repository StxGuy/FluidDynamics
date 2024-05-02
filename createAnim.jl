# Create a GIF animation file from the stillframes produced by the Fortran code
# Watch it using ]gifview -a test.gif

using Plots

anim = @animate for it in 1:9999

    file = open("video/test"*string(it)*".txt")
    lines = readlines(file)
    close(file)
    
    M = [parse.(Float64,split(line)) for line in lines]
    M = hcat(M...)
    heatmap(M)#,clim=(0,0.14))    
end

gif(anim,"test.gif",fps=500)
