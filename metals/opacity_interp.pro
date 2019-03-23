PRO opacity_interp

;Back-grids enriched C & O opacities to EZ format from new STARS format

CD, Current = currentDir

opacities_in=fltarr(31)
rhos_in=fltarr(31)
Rs = findgen(31)*0.5 - 7.0
temps = findgen(127)*0.05 + 3.0
rhos = findgen(90)*0.25 - 12.0
Xs = [0.0,0.03,0.1,0.35,0.7]
tempin = 1.0
opacin1 = 1.0
opacin2 = opacin1
opacin3 = opacin1
opacin4 = opacin1
opacin5 = opacin1
opacin6 = opacin1
opacin7 = opacin1
opacin8 = opacin1
opacin9 = opacin1
opacin10 = opacin1
opacin11 = opacin1
opacin12 = opacin1
opacin13 = opacin1
opacin14 = opacin1
opacin15 = opacin1
opacin16 = opacin1
opacin17 = opacin1
opacin18 = opacin1
opacin19 = opacin1
opacin20 = opacin1
opacin21 = opacin1
opacin22 = opacin1
opacin23 = opacin1
opacin24 = opacin1
opacin25 = opacin1
opacin26 = opacin1
opacin27 = opacin1
opacin28 = opacin1
opacin29 = opacin1
opacin30 = opacin1
opacin31 = opacin1
Xin = opacin1

OpenR,lun,'COtables_z000_noenrich.dat',/Get_Lun
OpenW,lun2,'EZ_EOS_prelim.dat',/Get_Lun

FOR i=0,4 DO BEGIN
    ReadF,lun,Xin,format='(F1.2)'
    IF Xin NE Xs(i) print 'Error 1'
    FOR j=0,90 DO BEGIN
       ReadF,lun,tempin,opacin1,opacin2,opacin3,opacin4,opacin5,
        $opacin6,opacin7,opacin8,opacin9,opacin10,
        $opacin11,opacin12,opacin13,opacin14,opacin15,
        $opacin16,opacin17,opacin18,opacin19,opacin20,
        $opacin21,opacin22,opacin23,opacin24,opacin25,
        $opacin26,opacin27,opacin28,opacin29,opacin30,opacin31,format='(F5.2, 31F7.3)
       IF tempin NE temps(i) print 'Error 2'
       opacities_in(0)=opacin1
       opacities_in(1)=opacin2
       opacities_in(2)=opacin3
       opacities_in(3)=opacin4
       opacities_in(4)=opacin5
       opacities_in(5)=opacin6
       opacities_in(6)=opacin7
       opacities_in(7)=opacin8
       opacities_in(8)=opacin9
       opacities_in(9)=opacin10
       opacities_in(10)=opacin11
       opacities_in(11)=opacin12
       opacities_in(12)=opacin13
       opacities_in(13)=opacin14
       opacities_in(14)=opacin15
       opacities_in(15)=opacin16
       opacities_in(16)=opacin17
       opacities_in(17)=opacin18
       opacities_in(18)=opacin19
       opacities_in(19)=opacin20
       opacities_in(20)=opacin21
       opacities_in(21)=opacin22
       opacities_in(22)=opacin23
       opacities_in(23)=opacin24
       opacities_in(24)=opacin25
       opacities_in(25)=opacin26
       opacities_in(26)=opacin27
       opacities_in(27)=opacin28
       opacities_in(28)=opacin29
       opacities_in(29)=opacin30
       opacities_in(30)=opacin31
       rhos_in = 

 lambda(i) = lambdat
    expot(i) = expott
    lgf(i) = lgft
    wobs(i) = wobst/10.
    abun(i) = abunt
    nor(i) = normt
    shiff(i) = shiftt
    ENDFOR
ENDFOR

Free_Lun, lun

IF KEYWORD_SET(mK) THEN wobs = 1.0e3/(1.0/lambda - 0.5e-10*wobs) - $
       1.0e3/(1.0/lambda + 0.5e-10*wobs)

IF plotcode EQ 1 THEN BEGIN
    IF NOT KEYWORD_SET(OPlot) THEN BEGIN
       IF KEYWORD_SET(xrange) THEN Plot, wobs, abun+baseabun1, /YNOZERO, YRANGE=yvals, $
          PSYM=8, YSTYLE=1, YTITLE = 'log !4e!X'+whatIsIt, XTITLE='Equivalent Width [pm]',$
          xrange=xrange, /nodata ELSE $
          Plot, wobs, abun+baseabun1, /YNOZERO, YRANGE=yvals,$
          PSYM=8, YSTYLE=1, YTITLE = 'log !4e!X'+whatIsIt, XTITLE='Equivalent Width [pm]', /nodata
    ENDIF
    OPlot, wobs, abun+baseabun1, PSYM=8, color=colour[0]
    fit = LinFit(wobs, abun+baseabun1)
ENDIF
IF plotcode EQ 2 THEN BEGIN
    IF NOT KEYWORD_SET(OPlot) THEN BEGIN
       IF KEYWORD_SET(xrange) THEN Plot, expot, abun+baseabun1, /YNOZERO, YRANGE=yvals, $
          PSYM=8, YSTYLE=1, YTITLE = 'log !4e!X'+whatIsIt, XTITLE='Lower Excitation Potential [eV]',$
          xrange=xrange, /nodata ELSE $
          Plot, expot, abun+baseabun1, /YNOZERO, YRANGE=yvals,$
          PSYM=8, YSTYLE=1, YTITLE = 'log !4e!X'+whatIsIt, XTITLE='Lower Excitation Potential [eV]', /nodata
    ENDIF
    OPlot, expot, abun+baseabun1, PSYM=8, color=colour[0]
    fit = LinFit(expot, abun+baseabun1)
ENDIF
IF plotcode EQ 3 THEN BEGIN
    IF NOT KEYWORD_SET(OPlot) THEN BEGIN
       IF KEYWORD_SET(xrange) THEN Plot, lambda, abun+baseabun1, /YNOZERO, YRANGE=yvals, $
          PSYM=8, YSTYLE=1, YTITLE = 'log !4e!X'+whatIsIt, XTITLE='Wavelength [nm]',$
          xrange=xrange, /nodata ELSE $
          Plot, lambda, abun+baseabun1, /YNOZERO, YRANGE=yvals, $
          PSYM=8, YSTYLE=1, YTITLE = 'log !4e!X'+whatIsIt, XTITLE='Wavelength [nm]', /nodata
    ENDIF
    OPlot, lambda, abun+baseabun1, PSYM=8, color=colour[0]
    fit = LinFit(lambda, abun+baseabun1)
ENDIF
IF NOT KEYWORD_SET(nofit) THEN Plots, [[!x.crange(0), fit(1)*!x.crange(0) + fit(0)], $
   [!x.crange(1), fit(1)*!x.crange(1) + fit(0)]], color=colour[0]

IF KEYWORD_SET(infile2) THEN BEGIN

    plotsym, points[1,0], points[1,1], fill=points[1,2]

    OpenR,lun,currentDir+'/'+infile2,/Get_Lun

    IF KEYWORD_SET(nlines2) THEN BEGIN
       nlines = nlines2
       lambda=DblArr(nlines)
       expot=FltArr(nlines)
       lgf=expot
       wobs=lambda
       abun=expot
       nor=expot
       shiff=expot
    ENDIF

    FOR i=0,nlines-1 DO BEGIN
        ReadF,lun,lambdat,expott,lgft,wobst,abunt,normt,shiftt,format='(F11.5,6F8.3)
       lambda(i) = lambdat
       expot(i) = expott
       lgf(i) = lgft
       wobs(i) = wobst/10.
       abun(i) = abunt
       nor(i) = normt
       shiff(i) = shiftt
    ENDFOR

    Free_Lun, lun

    IF KEYWORD_SET(mK) THEN wobs = 1.0e3/(1.0/lambda - 0.5e-10*wobs) - $
       1.0e3/(1.0/lambda + 0.5e-10*wobs)

    IF plotcode EQ 1 THEN BEGIN
       OPlot, wobs, abun+baseabun2, PSYM=8, color=colour[1]
       fit = LinFit(wobs, abun+baseabun2)
    ENDIF
    IF plotcode EQ 2 THEN BEGIN
       OPlot, expot, abun+baseabun2, PSYM=8, color=colour[1]
       fit = LinFit(expot, abun+baseabun2)
    ENDIF
    IF plotcode EQ 3 THEN BEGIN
       OPlot, lambda, abun+baseabun2, PSYM=8, color=colour[1]
       fit = LinFit(lambda, abun+baseabun2)
    ENDIF
    IF NOT KEYWORD_SET(nofit) THEN Plots, [[!x.crange(0), fit(1)*!x.crange(0) + fit(0)], $
       [!x.crange(1), fit(1)*!x.crange(1) + fit(0)]], LINESTYLE=2*$
       (1-((KEYWORD_SET(ps) AND colour[0] NE 0 AND colour[1] NE 0) OR $
           (NOT KEYWORD_SET(ps) AND colour[0] NE -1 AND colour[1] NE -1))), color=colour[1]

    IF KEYWORD_SET(model1) AND KEYWORD_SET(model2) THEN BEGIN
       models = [model1, model2]
       locateLegend, legendLocation, points, models, colour=colour
    ENDIF

ENDIF

