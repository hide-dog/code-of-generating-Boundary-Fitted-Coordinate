using Printf

println("\n")
println("------------------------------------------")
println("--  重複しないように                     --")
println("--  データセットを用意すること           --")
println("------------------------------------------")
println("\n")


function main()
    # 座標点数 ( 格子数 = nr-1 )
    nr    = 200
    nwing = 200
    
    
    fout     = "xy_NACA0012"         # output file name
    NACAname = "NACA0012"            # NACA0012 type
    r    = 2.5                       # O型格子の半径
    cx   = 0.5                       # 翼データはx=0~1のため、その間の0.5をO型の中心とする。
    dnum = 4                         # NACA0012 type = 4, NACA64A010 type = 6 
    
    
    #=
    fout     = "xy_hayabusa"         # output file name
    NACAname = "hayabusa_shape"      # NACA0012 type
    r    = 0.9                       # O型格子の半径
    cx   = -0.101                    # O型格子の中心
    dnum = 4                         # NACA0012 type = 4, NACA64A010 type = 6 
    =#
    
    #=
    fout     = "xy_NACA64A010"       # output file name
    NACAname = "NACA64A010"          # NACA64A010 type
    r    = 2.5                       # O型格子の半径
    cx   = 0.5                       # 翼データはx=0~1のため、その間の0.5をO型の中心とする。
    dnum = 6                         # NACA0012 type = 4, NACA64A010 type = 6 
    =#
    
    #=
    fout     = "xy_NACA4412"       # output file name
    NACAname = "NACA4412"          # NACA64A010 type
    r    = 2.5                       # O型格子の半径
    cx   = 0.5                       # 翼データはx=0~1のため、その間の0.5をO型の中心とする。
    dnum = 6                         # NACA0012 type = 4, NACA64A010 type = 6 
    =#
    
    
    println("\n")
    println("------------------------------------------")
    println("--  格子の切れ目は下記のデータを使用する --")
    println("--  4d : 下面の最後のデータ             --")
    println("--  6d : 一番最初のデータ               --")
    println("------------------------------------------")
    println("\n")
    
    if dnum == 4
        surface_len, edgex, edgey, aroundx, aroundy = read_NACA_4d(NACAname)
    elseif dnum == 6
        surface_len, edgex, edgey, aroundx, aroundy = read_NACA_6d(NACAname)
    else
        println("\n")
        println("--  データ読み込みの方法を決めてください --")
        println("--  datanum = 4 or 6                   --")
        println("\n")
        throw(UndefVarError(:x))
    end
    
    x1, y1 = cal_area1(edgex, edgey, r, cx, nr)
    x2, y2 = cal_area2(surface_len, aroundx, aroundy,nwing)
    x3, y3 = cal_area1(edgex, edgey, r, cx, nr)
    x4, y4 = cal_area4(edgex, edgey, r, cx, nwing)
    
    x,y = algebraic_grid_gene(nr,nwing,x1,x2,x3,x4,y1,y2,y3,y4)
    
    open(fout,"w") do f
        a1 = @sprintf("%8.8e", nr)
        a2 = @sprintf("%8.8e", nwing)
        write(f, a1*" "*a2*"\n")
        for i in 1:nr
            for j in 1:nwing
                ai = string(i)
                aj = string(j)
                a1 = @sprintf("%8.8e", x[i,j])
                a2 = @sprintf("%8.8e", y[i,j])
                write(f, ai*" "*aj*" "*a1*" "*a2*"\n")
            end
        end
    end

    open("External shape","w") do f
        for j in 1:nwing
            a1 = @sprintf("%8.8e", x2[j])
            a2 = @sprintf("%8.8e", y2[j])
            write(f, a1*" "*a2*"\n")
        end
    end
end

function cal_area1(edgex, edgey, r, cx, nr)
    x1 = zeros(nr)
    y1 = zeros(nr)
    
    dx  = (r - edgex + cx) / (nr-1)
    
    for i in 1:nr
        x1[i] = edgex + (i-1)*dx
        y1[i] = edgey
    end
    
    return x1,y1
end

function cal_area2(surface_len, aroundx, aroundy, nwin)
    x2 = zeros(nwin)
    y2 = zeros(nwin)

    dx  = surface_len/(nwin-1)
    naround = length(aroundx)

    pointlen  = 1e-50
    aroundlen = 0.0
    l = 0

    x2[1] = aroundx[1]
    y2[1] = aroundy[1]
    for i in 2:nwin
        pointlen  += dx
        while (pointlen - aroundlen)/(aroundlen) >= 1e-5
            l += 1
            aroundlen += ((aroundx[l]-aroundx[l+1])^2 + (aroundy[l]-aroundy[l+1])^2)^0.5
        end

        p1x = aroundx[l]
        p1y = aroundy[l]
        p2x = aroundx[l+1]
        p2y = aroundy[l+1]
        dis_2top = aroundlen - pointlen
        dis_2to1 = ((p1x-p2x)^2 + (p1y-p2y)^2)^0.5

        x2[i] = dis_2top/dis_2to1*p1x + (dis_2to1-dis_2top)/dis_2to1*p2x
        y2[i] = dis_2top/dis_2to1*p1y + (dis_2to1-dis_2top)/dis_2to1*p2y
    end
    return x2,y2
end

function cal_area4(edgex, edgey, r, cx, nwin)
    x4 = zeros(nwin)
    y4 = zeros(nwin)

    len = 2*pi
    dt  = len/(nwin-1)
    
    for i in 1:nwin
        t     = 2*pi - (i-1)*dt
        x4[i] = r * cos(t) + cx
        y4[i] = r * sin(t)
    end
    return x4,y4
end

# ----------------------------------------------
# -- algebraic_grid_gene(nx,ny,x1,x2,x3,x4,y1,y2,y3,y4)
# -- 幾何学的計算により、境界適合格子を作成
# ----------------------------------------------
function algebraic_grid_gene(nx,ny,x1,x2,x3,x4,y1,y2,y3,y4)
    x = zeros(nx,ny)
    y = zeros(nx,ny)
    
    for j in 2:nx-1
        for i in 2:ny-1
            xi  = (i-1)/(ny-1)
            eta = (j-1)/(nx-1)
            x[i,j] = (1-eta)*x1[i] + eta*x3[i] + (1-xi)*x2[j] + xi*x4[j] - 
                    (xi*eta*x3[nx] + xi*(1-eta)*x1[nx] + eta*(1-xi)*x3[1] + (1-xi)*(1-eta)*x1[1])
            y[i,j] = (1-eta)*y1[i] + eta*y3[i] + (1-xi)*y2[j] + xi*y4[j] - 
                    (xi*eta*y3[ny] + xi*(1-eta)*y1[ny] + eta*(1-xi)*y3[1] + (1-xi)*(1-eta)*y1[1])
        end
    end
    
    # Substitution of bd
    for i in 1:nx
        x[i,1]  = x1[i]
        y[i,1]  = y1[i]
        x[i,ny] = x3[i]
        y[i,ny] = y3[i]
    end
    for j in 1:ny
        x[1,j]  = x2[j]
        y[1,j]  = y2[j]
        x[nx,j] = x4[j]
        y[nx,j] = y4[j]
    end

    return x,y
end

# ----------------------------------------------
# -- read_NACA_4d(NACAname)
# -- データ読み込み、(NACA0012のデータセットの読み方)
# ----------------------------------------------
function read_NACA_4d(NACAname)
    
    fff=[]
    open(NACAname, "r") do f                  # 全行格納
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)     # 改行分割(\n)
    
    for i in 1:length(fff)
        fff[i] = replace(fff[i]," \r" => "")  # 改行削除(\r)
    end
    
    # 上面、下面データ数の読み込み
    temp   = split(fff[2]," ")
    filter!(e->e≠"",temp)

    nt = Int(parse(Float64,temp[1]))
    nb = Int(parse(Float64,temp[2]))
    
    # 上面、下面ポイントの読み込み
    top = zeros(nt, 2)
    bot = zeros(nb, 2)

    skip = 2
    for i in 1:nt
        temp = split(fff[i+skip]," ")
        filter!(e->e≠"",temp)
        
        top[i,1] = parse(Float64,temp[1])
        top[i,2] = parse(Float64,temp[2])
    end

    skip = 2 + nt
    for i in 1:nb
        temp = split(fff[i+skip]," ")
        filter!(e->e≠"",temp)

        bot[i,1] = parse(Float64,temp[1])
        bot[i,2] = parse(Float64,temp[2])
    end

    # 表面データの結合
    aroundx = zeros(nt+nb+1)
    aroundy = zeros(nt+nb+1)
    for i in 1:nb
        aroundx[i] = bot[nb+1-i,1]
        aroundy[i] = bot[nb+1-i,2]
    end
    for i in 1:nt
        aroundx[i+nb] = top[i,1]
        aroundy[i+nb] = top[i,2]
    end
    aroundx[nt+nb+1] = bot[nb,1]
    aroundy[nt+nb+1] = bot[nb,2]
    
    # 表面線長さ、エッジ座標の登録
    surface_len = 0.0
    edgex = 0.0
    edgey = 0.0
    
    n = nt+nb
    temp = 0.0
    for i in 1:n
        temp += ((aroundx[i]-aroundx[i+1])^2 + (aroundy[i]-aroundy[i+1])^2)^0.5
    end
    surface_len = temp

    edgex = aroundx[1]
    edgey = aroundy[1]
        
    return surface_len, edgex, edgey, aroundx, aroundy
end

# ----------------------------------------------
# -- read_NACA_6d(NACAname)
# -- データ読み込み、(NACA64A010のデータセットの読み方)
# ----------------------------------------------
function read_NACA_6d(NACAname)
    
    fff=[]
    open(NACAname, "r") do f                  # 全行格納
        fff = read(f,String)
    end 
    fff = split(fff,"\n",keepempty=false)     # 改行分割(\n)
    
    for i in 1:length(fff)
        fff[i] = replace(fff[i]," \r" => "")  # 改行削除(\r)
    end
    
    # 表面データの読み込み
    skip = 1
    n    = length(fff)-skip
    aroundx = zeros(n+1)
    aroundy = zeros(n+1)
    taroundx = zeros(n+1)
    taroundy = zeros(n+1)
    for i in 1:n
        temp = split(fff[i+skip]," ")
        filter!(e->e≠"",temp)
        
        taroundx[i] = parse(Float64,temp[1])
        taroundy[i] = parse(Float64,temp[2])
    end
    taroundx[n+1] = taroundx[1]
    taroundy[n+1] = taroundy[1]

    for i in 1:n+1
        aroundx[i] = taroundx[n+1-(i-1)]
        aroundy[i] = taroundy[n+1-(i-1)]
    end

    # 表面線長さ、エッジ座標の登録
    surface_len = 0.0
    edgex = 0.0
    edgey = 0.0

    temp = 0.0
    for i in 1:n
        temp += ((aroundx[i]-aroundx[i+1])^2 + (aroundy[i]-aroundy[i+1])^2)^0.5
    end
    surface_len = temp

    edgex = aroundx[1]
    edgey = aroundy[1]
    
    return surface_len, edgex, edgey, aroundx, aroundy
end

# ---------------------------------
main() 