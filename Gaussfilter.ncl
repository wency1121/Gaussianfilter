
; M-term Gaussian-Type Filter
; Input variables: x(n), m
; m: the term number used to running mean
; it must be an odd number.
; Output variables: y(n)
; y: the filtered series of x.
; Work parameters and array: c, cgm and ck(-(m-1)/2:(m-1)/2)
; c: a tunable parameter, generally, c>2.0.
; cgm: variance of Gaussian distribution.
; ck: weight


undef("gaussfiltern")
function gaussfiltern(x,m)
begin
    ; Get the input array's dimension sizes and rank.
    x_dims = dimsizes(x)
    x_rank = dimsizes(x_dims)

    ; Get the number of time steps.
    n = x_dims(0)
    dims=n+m-1

    pi= 3.1415926
    c=2.15
    nl=(m-1)/2
    cgm=tofloat(nl)/c
    c1=1./(cgm*sqrt(2.*pi))
    
    wgt=new((/m/),"float")
    wgt(nl)=c1
    a=ispan(1,nl,1)
    wgt(nl+1:2*nl)=c1*exp(-a*a/(2.*cgm*cgm))
    wgt(0:nl-1:-1)=wgt(nl+1:2*nl)
    
    if(x_rank.eq.1)then
        xw = new(dims, "float", getFillValue(x))
        xw(0:nl)=x(0)
        xw(nl+1:n+nl-2)=x(1:n-2)
        xw(n+nl-1:dims-1)=x(n-1)
        y=wgt_runave_n(xw, wgt, 0, 0)
        y:=y(nl:dims-nl-1)
    end if

    if(x_rank.eq.2)then
        dims1 = x_dims(1)
        xw    = new((/dims, dims1/), "float", getFillValue(x))
        xw(0:nl,:)          = conform(xw(0:nl,:),x(0,:),(/1/))
        xw(nl+1:n+nl-2,:)   = x(1:n-2,:)
        xw(n+nl-1:dims-1,:) = conform(xw(n+nl-1:dims-1,:),x(n-1,:),(/1/))
        y  = wgt_runave_n(xw, wgt, 0, 0)
        y := y(nl:dims-nl-1,:)
    end if

    if(x_rank.eq.3)then
        dims1 = x_dims(1)
        dims2 = x_dims(2)
        xw    = new((/dims, dims1, dims2/), "float", getFillValue(x))
        xw(0:nl,:,:)          = conform(xw(0:nl,:,:),x(0,:,:),(/1,2/))
        xw(nl+1:n+nl-2,:,:)   = x(1:n-2,:,:)
        xw(n+nl-1:dims-1,:,:) = conform(xw(n+nl-1:dims-1,:,:),x(n-1,:,:),(/1,2/))
        y  = wgt_runave_n(xw, wgt, 0, 0)
        y := y(nl:dims-nl-1,:,:)
    end if

    if(x_rank.eq.4)then
        dims1 = x_dims(1)
        dims2 = x_dims(2)
        dims3 = x_dims(3)
        xw    = new((/dims, dims1, dims2, dims3/), "float", getFillValue(x))
        xw(0:nl,:,:,:)          = conform(xw(0:nl,:,:,:),x(0,:,:,:),(/1,2,3/))
        xw(nl+1:n+nl-2,:,:,:)   = x(1:n-2,:,:,:)
        xw(n+nl-1:dims-1,:,:,:) = conform(xw(n+nl-1:dims-1,:,:,:),x(n-1,:,:,:),(/1,2,3/))
        y  = wgt_runave_n(xw, wgt, 0, 0)
        y := y(nl:dims-nl-1,:,:,:)
    end if

    if(x_rank.eq.5)then
        dims1 = x_dims(1)
        dims2 = x_dims(2)
        dims3 = x_dims(3)
        dims4 = x_dims(4)
        xw    = new((/dims, dims1, dims2, dims3, dims4/), "float", getFillValue(x))
        xw(0:nl,:,:,:,:)          = conform(xw(0:nl,:,:,:,:),x(0,:,:,:,:),(/1,2,3,4/))
        xw(nl+1:n+nl-2,:,:,:,:)   = x(1:n-2,:,:,:,:)
        xw(n+nl-1:dims-1,:,:,:,:) = conform(xw(n+nl-1:dims-1,:,:,:,:),x(n-1,:,:,:,:),(/1,2,3,4/))
        y  = wgt_runave_n(xw, wgt, 0, 0)
        y := y(nl:dims-nl-1,:,:,:,:)
    end if

    copy_VarCoords(x, y)

    return(y)
end

