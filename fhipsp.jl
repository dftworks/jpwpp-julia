function read_in_psp_fhi()
    lines = readlines("14-Si.LDA.fhi")
    zatom, zion = parse.(Float64, split(lines[2])[1:2])
    pspcod, pspxc, lmax, lloc, mmax = parse.(Int32, split(lines[3])[1:5])
    
end

read_in_psp_fhi()
