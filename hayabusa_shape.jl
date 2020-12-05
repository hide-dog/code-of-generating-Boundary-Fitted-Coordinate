# ----------------------------------
# -- Reference
# Planetary Mission Entry Vehicles Quick Reference Guide. Version 3.0
# Davies, Carol (Eloret Corp. Moffett Field, CA, United States), Arcadi, Marla
# January 1, 2006
# https://ntrs.nasa.gov/citations/20070022789
# ----------------------------------

using Printf

function cal_area1()
    n = 20
    r = 0.202
    finx, finy = -0.15724, 0.12971
    
    x1 = zeros(n)
    y1 = zeros(n)

    t = asin(finy/r)
    dt = t/n

    for i in 1:n
        x1[i] = -r * cos(dt*(i-1))
        y1[i] =  r * sin(dt*(i-1))
    end
    return x1,y1
end

function cal_area2()
    n = 5
    stax, stay = -0.15724, 0.12971
    finx, finy = -0.085, 0.202
    
    x2 = zeros(n)
    y2 = zeros(n)
    
    l   = ((stax-finx)^2 + (stay-finy)^2)^0.5
    dx  = l/n
    
    templ = 0.0
    for i in 1:n
        x2[i] = (l-templ)/l*stax + templ/l*finx
        y2[i] = (l-templ)/l*stay + templ/l*finy
        templ += dx
    end
    return x2,y2
end

function cal_area3()
    n = 5
    stax, stay = -0.085, 0.202
    finx, finy = 0.0, 0.117
    
    x3 = zeros(n)
    y3 = zeros(n)
    
    l   = ((stax-finx)^2 + (stay-finy)^2)^0.5
    dx  = l/n
    
    templ = 0.0
    for i in 1:n
        x3[i] = (l-templ)/l*stax + templ/l*finx
        y3[i] = (l-templ)/l*stay + templ/l*finy
        templ += dx
    end
    return x3,y3
end

function cal_area4()
    n = 5
    stax, stay = 0.0, 0.117
    finx, finy = 0.0, 0.0
    
    x4 = zeros(n+1)
    y4 = zeros(n+1)
    
    l   = ((stax-finx)^2 + (stay-finy)^2)^0.5
    dx  = l/n
    
    templ = 0.0
    for i in 1:n
        x4[i] = (l-templ)/l*stax + templ/l*finx
        y4[i] = (l-templ)/l*stay + templ/l*finy
        templ += dx
    end
    x4[n+1] = finx
    y4[n+1] = finy
    
    return x4,y4
end

function matome(x1,x2,x3,x4,y1,y2,y3,y4)
    n = length(x1) + length(x2) + length(x3) + length(x4)

    x = zeros(n)
    y = zeros(n)

    n1 = length(x1)
    n2 = length(x1) + length(x2)
    n3 = length(x1) + length(x2) + length(x3)
    n4 = length(x1) + length(x2) + length(x3) + length(x4)
    for i in 1:n1
        x[i] = x1[i]
        y[i] = y1[i]
    end
    for i in n1+1:n2
        x[i] = x2[i-n1]
        y[i] = y2[i-n1]
    end
    for i in n2+1:n3
        x[i] = x3[i-n2]
        y[i] = y3[i-n2]
    end
    for i in n3+1:n4
        x[i] = x4[i-n3]
        y[i] = y4[i-n3]
    end
    return x, y     
end

function cal_bot(topx, topy)
    n = length(topx)
    botx = zeros(n)
    boty = zeros(n)

    for i in 1:n
        botx[i] = topx[i]
        boty[i] = -topy[i]
    end

    tx = zeros(n-2)
    ty = zeros(n-2)
    for i in 1:n-2
        tx[i] = topx[i+1]
        ty[i] = topy[i+1]
    end
    
    return tx,ty,botx,boty
end

function main()
    fout = "hayabusa_shape"
        
    x1, y1 = cal_area1()
    x2, y2 = cal_area2()
    x3, y3 = cal_area3()
    x4, y4 = cal_area4()
    
    topx, topy = matome(x1,x2,x3,x4,y1,y2,y3,y4)
    topx, topy, botx, boty = cal_bot(topx, topy)
    
    nt = length(topx)
    nb = length(botx)
    open(fout,"w") do f
        write(f,"hayabusa \n")
        a1 = string(nt)
        a2 = string(nb)
        write(f, a1*" "*a2*"\n")
        write(f, "\n")

        for i in 1:nt
            a1 = @sprintf("%8.8e", topx[i])
            a2 = @sprintf("%8.8e", topy[i])
            write(f, a1*" "*a2*"\n")
        end
        write(f, "\n")
        for i in 1:nb
            a1 = @sprintf("%8.8e", botx[i])
            a2 = @sprintf("%8.8e", boty[i])
            write(f, a1*" "*a2*"\n")
        end
    end
end

# ---------------------------------
main() 