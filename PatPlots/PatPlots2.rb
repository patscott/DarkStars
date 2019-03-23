class PatPlots2

    include Math
    include Tioga
    include FigureConstants
       
    Metallicity = 0.01
    Epsilon = 1e-6
    HREpsilon = 0.1
              attr_accessor :metal, :mass, :rhowimp, :v_star, :t1, :tcmin, :lwimp0, :lwimp1, :lpp0, :lpp1, :lCNO0, :lCNO1, :cap, :cond_eff, :ms_Age, :conv_M_s1, :conv_M_e1, :conv_M_s2, :conv_M_e2, :conv_M_s3, :conv_M_e3
attr_accessor :mchi_2, :sigmav_2, :mass_2, :rhowimp_2, :v_star_2, :t1_2, :tcmin_2, :lwimp0_2, :lwimp1_2, :lpp0_2, :lpp1_2, :lCNO0_2, :lCNO1_2, :cap_2, :cond_eff_2, :ms_Age_2, :conv_M_s1_2, :conv_M_e1_2, :conv_M_s2_2, :conv_M_e2_2, :conv_M_s3_2, :conv_M_e3_2
attr_accessor :mchi_3, :sigmav_3, :mass_3, :rhowimp_3, :v_star_3, :t1_3, :tcmin_3, :lwimp0_3, :lwimp1_3, :lpp0_3, :lpp1_3, :lCNO0_3, :lCNO1_3, :cap_3, :cond_eff_3, :ms_Age_3, :conv_M_s1_3, :conv_M_e1_3, :conv_M_s2_3, :conv_M_e2_3, :conv_M_s3_3, :conv_M_e3_3
attr_accessor :metal_4, :mass_4, :rhowimp_4, :v_star_4, :t1_4, :tcmin_4, :lwimp0_4, :lwimp1_4, :lpp0_4, :lpp1_4, :lCNO0_4, :lCNO1_4, :cap_4, :cond_eff_4, :ms_Age_4, :conv_M_s1_4, :conv_M_e1_4, :conv_M_s2_4, :conv_M_e2_4, :conv_M_s3_4, :conv_M_e3_4
attr_accessor :hrdata
attr_accessor :galr_6, :mass_6, :rhowimp_6, :v_star_6, :t1_6, :tcmin_6, :lwimp0_6, :lwimp1_6, :lpp0_6, :lpp1_6, :lCNO0_6, :lCNO1_6, :cap_6, :cond_eff_6, :ms_Age_6, :conv_M_s1_6, :conv_M_e1_6, :conv_M_s2_6, :conv_M_e2_6, :conv_M_s3_6, :conv_M_e3_6
attr_accessor :galr_7, :mass_7, :rhowimp_7, :v_star_7, :t1_7, :tcmin_7, :lwimp0_7, :lwimp1_7, :lpp0_7, :lpp1_7, :lCNO0_7, :lCNO1_7, :cap_7, :cond_eff_7, :ms_Age_7, :conv_M_s1_7, :conv_M_e1_7, :conv_M_s2_7, :conv_M_e2_7, :conv_M_s3_7, :conv_M_e3_7
attr_accessor :galr_8, :mass_8, :rhowimp_8, :v_star_8, :t1_8, :tcmin_8, :lwimp0_8, :lwimp1_8, :lpp0_8, :lpp1_8, :lCNO0_8, :lCNO1_8, :cap_8, :cond_eff_8, :ms_Age_8, :conv_M_s1_8, :conv_M_e1_8, :conv_M_s2_8, :conv_M_e2_8, :conv_M_s3_8, :conv_M_e3_8
attr_accessor :orbit_9, :e_9, :mass_9, :rhowimp_9, :v_star_9, :t1_9, :tcmin_9, :lwimp0_9, :lwimp1_9, :lpp0_9, :lpp1_9, :lCNO0_9, :lCNO1_9, :cap_9, :cond_eff_9, :ms_Age_9, :conv_M_s1_9, :conv_M_e1_9, :conv_M_s2_9, :conv_M_e2_9, :conv_M_s3_9, :conv_M_e3_9
attr_accessor :orbit_10, :e_10, :mass_10, :rhowimp_10, :v_star_10, :t1_10, :tcmin_10, :lwimp0_10, :lwimp1_10, :lpp0_10, :lpp1_10, :lCNO0_10, :lCNO1_10, :cap_10, :cond_eff_10, :ms_Age_10, :conv_M_s1_10, :conv_M_e1_10, :conv_M_s2_10, :conv_M_e2_10, :conv_M_s3_10, :conv_M_e3_10
attr_accessor :orbit_11, :e_11, :mass_11, :rhowimp_11, :v_star_11, :t1_11, :tcmin_11, :lwimp0_11, :lwimp1_11, :lpp0_11, :lpp1_11, :lCNO0_11, :lCNO1_11, :cap_11, :cond_eff_11, :ms_Age_11, :conv_M_s1_11, :conv_M_e1_11, :conv_M_s2_11, :conv_M_e2_11, :conv_M_s3_11, :conv_M_e3_11
attr_accessor :orbit_12, :e_12, :mass_12, :rhowimp_12, :v_star_12, :t1_12, :tcmin_12, :lwimp0_12, :lwimp1_12, :lpp0_12, :lpp1_12, :lCNO0_12, :lCNO1_12, :cap_12, :cond_eff_12, :ms_Age_12, :conv_M_s1_12, :conv_M_e1_12, :conv_M_s2_12, :conv_M_e2_12, :conv_M_s3_12, :conv_M_e3_12
attr_accessor :orbit_13, :e_13, :mass_13, :rhowimp_13, :v_star_13, :t1_13, :tcmin_13, :lwimp0_13, :lwimp1_13, :lpp0_13, :lpp1_13, :lCNO0_13, :lCNO1_13, :cap_13, :cond_eff_13, :ms_Age_13, :conv_M_s1_13, :conv_M_e1_13, :conv_M_s2_13, :conv_M_e2_13, :conv_M_s3_13, :conv_M_e3_13
attr_accessor :etrans14, :radius14

    def initialize
        patPlots2_init
	t.tex_preamble = t.tex_preamble + "\n\t\\usepackage{amssymb}\n"
        return self
    end
    
    def t
        @figure_maker
    end

    def colours
        @colours
    end

    def linetypes_2
        @linetypes_2
    end

    def linetypes_3
        @linetypes_3
    end

    def orbits
        @orbits
    end

    def background
        t.fill_color = White
        t.fill_frame
    end

    def patPlots2_init
        @figure_maker = FigureMaker.default
        t.def_eval_function { |str| eval(str) }

        t.save_dir = 'PatPlots_out'
        t.def_figure("etrans"){plot_etrans}
        t.def_figure("CapPlot") {plot_cap}
        t.def_figure("LwimpPlot") {plot_Lwimp}
        t.def_figure("AnnPlot") {plot_ann}
        t.def_figure("MwimpPlot") {plot_Mwimp}
        t.def_figure("HRPlotp1") {plot_HRPlotp1}
        t.def_figure("TrhoPlotp1") {plot_Trhop1}
        t.def_figure("HRPlot0") {plot_HRPlot0}
        t.def_figure("TrhoPlot0") {plot_Trho0}
        t.def_figure("HRPlotn1") {plot_HRPlotn1}
        t.def_figure("TrhoPlotn1") {plot_Trhon1}
        t.def_figure("HRPlotn3") {plot_HRPlotn3}
        t.def_figure("TrhoPlotn3") {plot_Trhon3}
        t.def_figure("TcminPlot") {plot_Tcmin}
        t.def_figure("CondPlot") {plot_cond}
        t.def_figure("ppPlot") {plot_pp}
        t.def_figure("CNOPlot") {plot_CNO}
        t.def_figure("ConvPlot03") {plot_conv03}
        t.def_figure("ConvPlot04") {plot_conv04}
        t.def_figure("ConvPlot05") {plot_conv05}
        t.def_figure("ConvPlot06") {plot_conv06}
        t.def_figure("ConvPlot07") {plot_conv07}
        t.def_figure("ConvPlot08") {plot_conv08}
        t.def_figure("ConvPlot09") {plot_conv09}
        t.def_figure("ConvPlot10") {plot_conv10}
        t.def_figure("ConvPlot12") {plot_conv12}
        t.def_figure("ConvPlot14") {plot_conv14}
        t.def_figure("ConvPlot16") {plot_conv16}
        t.def_figure("ConvPlot18") {plot_conv18}
        t.def_figure("ConvPlot20") {plot_conv20}
        t.def_figure("TMSPlot") {plot_TMS}
        t.def_figure("circADPlot") {plot_circAD}
        t.def_figure("circNFWPlot") {plot_circNFW}
	t.def_figure("elliptical10yr") {plot_elliptical_10yr}
	t.def_figure("elliptical50yr") {plot_elliptical_50yr}
	t.def_figure("ellipticalcpc") {plot_elliptical_cpc}
	t.def_figure("ellipticalhalo") {plot_elliptical_halo}
	t.def_figure("ellipticallwimp") {plot_elliptical_lwimp}

        @data = [ # columns in ZAMSgridII.out
            @metal = Dvector.new,
            @mass = Dvector.new,
            @rhowimp = Dvector.new,
            @v_star = Dvector.new,
            @t1 = Dvector.new,
            @tcmin = Dvector.new,
            @lwimp0 = Dvector.new,
            @lwimp1 = Dvector.new,
            @lpp0 = Dvector.new,
            @lpp1 = Dvector.new,
            @lCNO0 = Dvector.new,
            @lCNO1 = Dvector.new,
            @cap = Dvector.new,
            @cond_eff = Dvector.new,
            @ms_Age = Dvector.new,
            @conv_M_s1 = Dvector.new,
            @conv_M_e1 = Dvector.new,
            @conv_M_s2 = Dvector.new,
            @conv_M_e2 = Dvector.new,
            @conv_M_s3 = Dvector.new,
            @conv_M_e3 = Dvector.new ]
        Dvector.read('ZAMSgridII.out', @data)     

        @data2 = [ # columns in m_chi.out
            @mchi_2 = Dvector.new,
            @sigmav_2 = Dvector.new,
            @mass_2 = Dvector.new,
            @rhowimp_2 = Dvector.new,
            @v_star_2 = Dvector.new,
            @t1_2 = Dvector.new,
            @tcmin_2 = Dvector.new,
            @lwimp0_2 = Dvector.new,
            @lwimp1_2 = Dvector.new,
            @lpp0_2 = Dvector.new,
            @lpp1_2 = Dvector.new,
            @lCNO0_2 = Dvector.new,
            @lCNO1_2 = Dvector.new,
            @cap_2 = Dvector.new,
            @cond_eff_2 = Dvector.new,
            @ms_Age_2 = Dvector.new,
            @conv_M_s1_2 = Dvector.new,
            @conv_M_e1_2 = Dvector.new,
            @conv_M_s2_2 = Dvector.new,
            @conv_M_e2_2 = Dvector.new,
            @conv_M_s3_2 = Dvector.new,
            @conv_M_e3_2 = Dvector.new ]
        Dvector.read('m_chi.out', @data2)     

        @data3 = [ # columns in sigmav.out
            @mchi_3 = Dvector.new,
            @sigmav_3 = Dvector.new,
            @mass_3 = Dvector.new,
            @rhowimp_3 = Dvector.new,
            @v_star_3 = Dvector.new,
            @t1_3 = Dvector.new,
            @tcmin_3 = Dvector.new,
            @lwimp0_3 = Dvector.new,
            @lwimp1_3 = Dvector.new,
            @lpp0_3 = Dvector.new,
            @lpp1_3 = Dvector.new,
            @lCNO0_3 = Dvector.new,
            @lCNO1_3 = Dvector.new,
            @cap_3 = Dvector.new,
            @cond_eff_3 = Dvector.new,
            @ms_Age_3 = Dvector.new,
            @conv_M_s1_3 = Dvector.new,
            @conv_M_e1_3 = Dvector.new,
            @conv_M_s2_3 = Dvector.new,
            @conv_M_e2_3 = Dvector.new,
            @conv_M_s3_3 = Dvector.new,
            @conv_M_e3_3 = Dvector.new ]
        Dvector.read('sigmav.out', @data3)     

        @data4 = [ # columns in HR.out
            @metal_4 = Dvector.new,
            @mass_4 = Dvector.new,
            @rhowimp_4 = Dvector.new,
            @v_star_4 = Dvector.new,
            @t1_4 = Dvector.new,
            @tcmin_4 = Dvector.new,
            @lwimp0_4 = Dvector.new,
            @lwimp1_4 = Dvector.new,
            @lpp0_4 = Dvector.new,
            @lpp1_4 = Dvector.new,
            @lCNO0_4 = Dvector.new,
            @lCNO1_4 = Dvector.new,
            @cap_4 = Dvector.new,
            @cond_eff_4 = Dvector.new,
            @ms_Age_4 = Dvector.new,
            @conv_M_s1_4 = Dvector.new,
            @conv_M_e1_4 = Dvector.new,
            @conv_M_s2_4 = Dvector.new,
            @conv_M_e2_4 = Dvector.new,
            @conv_M_s3_4 = Dvector.new,
            @conv_M_e3_4 = Dvector.new ]
        Dvector.read('HR.out', @data4)
	
	@data5 = [ # columns in z0*rho*M*.out
		@l_in = Dvector.new,
		@t_in = Dvector.new,
		@tc_in = Dvector.new,
		@rhoc_in = Dvector.new ]

	@points_in = Array.new

	@hrdata = Array.new
	@metal_4.each_index do |i|
		tempZ = @metal_4[i].to_s.sub("0.", "")
	        temprho = @rhowimp_4[i].to_s.sub(".0", "")
		tempmass = @mass_4[i].to_s.sub(".0", "")
	        Dvector.read('z'+tempZ+'rho'+temprho+'M'+tempmass+'ADMB.out', @data5)
                if (tempmass == "1.4") then
			@points_in = IO.readlines('z'+tempZ+'rho'+temprho+'M'+tempmass+'ADMB_points.out')
			@points_in.each_index do |j| @points_in[j] = @points_in[j][1,6].to_i end
                else @points_in = [1] end
		@hrdata[i] = Array.new
		@hrdata[i] = [Dvector.new(@l_in), Dvector.new(@t_in), Dvector.new(@tc_in), Dvector.new(@rhoc_in), @mass_4[i], Dvector.new(@points_in)]
	end

        @data6 = [ # columns in circorbitsAD.out
            @galr_6 = Dvector.new,
            @mass_6 = Dvector.new,
            @rhowimp_6 = Dvector.new,
            @v_star_6 = Dvector.new,
            @t1_6 = Dvector.new,
            @tcmin_6 = Dvector.new,
            @lwimp0_6 = Dvector.new,
            @lwimp1_6 = Dvector.new,
            @lpp0_6 = Dvector.new,
            @lpp1_6 = Dvector.new,
            @lCNO0_6 = Dvector.new,
            @lCNO1_6 = Dvector.new,
            @cap_6 = Dvector.new,
            @cond_eff_6 = Dvector.new,
            @ms_Age_6 = Dvector.new,
            @conv_M_s1_6 = Dvector.new,
            @conv_M_e1_6 = Dvector.new,
            @conv_M_s2_6 = Dvector.new,
            @conv_M_e2_6 = Dvector.new,
            @conv_M_s3_6 = Dvector.new,
            @conv_M_e3_6 = Dvector.new ]
        Dvector.read('circorbitsAD.out', @data6)

        @data7 = [ # columns in circorbitsNFW.out
            @galr_7 = Dvector.new,
            @mass_7 = Dvector.new,
            @rhowimp_7 = Dvector.new,
            @v_star_7 = Dvector.new,
            @t1_7 = Dvector.new,
            @tcmin_7 = Dvector.new,
            @lwimp0_7 = Dvector.new,
            @lwimp1_7 = Dvector.new,
            @lpp0_7 = Dvector.new,
            @lpp1_7 = Dvector.new,
            @lCNO0_7 = Dvector.new,
            @lCNO1_7 = Dvector.new,
            @cap_7 = Dvector.new,
            @cond_eff_7 = Dvector.new,
            @ms_Age_7 = Dvector.new,
            @conv_M_s1_7 = Dvector.new,
            @conv_M_e1_7 = Dvector.new,
            @conv_M_s2_7 = Dvector.new,
            @conv_M_e2_7 = Dvector.new,
            @conv_M_s3_7 = Dvector.new,
            @conv_M_e3_7 = Dvector.new ]
        Dvector.read('circorbitsNFW.out', @data7)

        @data8 = [ # columns in circorbitsAD_esc.out
            @galr_8 = Dvector.new,
            @mass_8 = Dvector.new,
            @rhowimp_8 = Dvector.new,
            @v_star_8 = Dvector.new,
            @t1_8 = Dvector.new,
            @tcmin_8 = Dvector.new,
            @lwimp0_8 = Dvector.new,
            @lwimp1_8 = Dvector.new,
            @lpp0_8 = Dvector.new,
            @lpp1_8 = Dvector.new,
            @lCNO0_8 = Dvector.new,
            @lCNO1_8 = Dvector.new,
            @cap_8 = Dvector.new,
            @cond_eff_8 = Dvector.new,
            @ms_Age_8 = Dvector.new,
            @conv_M_s1_8 = Dvector.new,
            @conv_M_e1_8 = Dvector.new,
            @conv_M_s2_8 = Dvector.new,
            @conv_M_e2_8 = Dvector.new,
            @conv_M_s3_8 = Dvector.new,
            @conv_M_e3_8 = Dvector.new ]
        Dvector.read('circorbitsAD_esc.out', @data8)

        @data9 = [ # columns in ellipticalAD.out
            @orbit_9 = Dvector.new,
	    @e_9 = Dvector.new,
            @mass_9 = Dvector.new,
            @rhowimp_9 = Dvector.new,
            @v_star_9 = Dvector.new,
            @t1_9 = Dvector.new,
            @tcmin_9 = Dvector.new,
            @lwimp0_9 = Dvector.new,
            @lwimp1_9 = Dvector.new,
            @lpp0_9 = Dvector.new,
            @lpp1_9 = Dvector.new,
            @lCNO0_9 = Dvector.new,
            @lCNO1_9 = Dvector.new,
            @cap_9 = Dvector.new,
            @cond_eff_9 = Dvector.new,
            @ms_Age_9 = Dvector.new,
            @conv_M_s1_9 = Dvector.new,
            @conv_M_e1_9 = Dvector.new,
            @conv_M_s2_9 = Dvector.new,
            @conv_M_e2_9 = Dvector.new,
            @conv_M_s3_9 = Dvector.new,
            @conv_M_e3_9 = Dvector.new ]
        Dvector.read('ellipticalAD.out', @data9)

        @data10 = [ # columns in ellipticalNFW.out
            @orbit_10 = Dvector.new,
	    @e_10 = Dvector.new,
            @mass_10 = Dvector.new,
            @rhowimp_10 = Dvector.new,
            @v_star_10 = Dvector.new,
            @t1_10 = Dvector.new,
            @tcmin_10 = Dvector.new,
            @lwimp0_10 = Dvector.new,
            @lwimp1_10 = Dvector.new,
            @lpp0_10 = Dvector.new,
            @lpp1_10 = Dvector.new,
            @lCNO0_10 = Dvector.new,
            @lCNO1_10 = Dvector.new,
            @cap_10 = Dvector.new,
            @cond_eff_10 = Dvector.new,
            @ms_Age_10 = Dvector.new,
            @conv_M_s1_10 = Dvector.new,
            @conv_M_e1_10 = Dvector.new,
            @conv_M_s2_10 = Dvector.new,
            @conv_M_e2_10 = Dvector.new,
            @conv_M_s3_10 = Dvector.new,
            @conv_M_e3_10 = Dvector.new ]
        Dvector.read('ellipticalNFW.out', @data10)

	@data11 = [ # columns in ellipticalAD_alt.out
            @orbit_11 = Dvector.new,
	    @e_11 = Dvector.new,
            @mass_11 = Dvector.new,
            @rhowimp_11 = Dvector.new,
            @v_star_11 = Dvector.new,
            @t1_11 = Dvector.new,
            @tcmin_11 = Dvector.new,
            @lwimp0_11 = Dvector.new,
            @lwimp1_11 = Dvector.new,
            @lpp0_11 = Dvector.new,
            @lpp1_11 = Dvector.new,
            @lCNO0_11 = Dvector.new,
            @lCNO1_11 = Dvector.new,
            @cap_11 = Dvector.new,
            @cond_eff_11 = Dvector.new,
            @ms_Age_11 = Dvector.new,
            @conv_M_s1_11 = Dvector.new,
            @conv_M_e1_11 = Dvector.new,
            @conv_M_s2_11 = Dvector.new,
            @conv_M_e2_11 = Dvector.new,
            @conv_M_s3_11 = Dvector.new,
            @conv_M_e3_11 = Dvector.new ]
        Dvector.read('ellipticalAD_alt.out', @data11)

	@data12 = [ # columns in ellipticalAD_alt2.out
            @orbit_12 = Dvector.new,
	    @e_12 = Dvector.new,
            @mass_12 = Dvector.new,
            @rhowimp_12 = Dvector.new,
            @v_star_12 = Dvector.new,
            @t1_12 = Dvector.new,
            @tcmin_12 = Dvector.new,
            @lwimp0_12 = Dvector.new,
            @lwimp1_12 = Dvector.new,
            @lpp0_12 = Dvector.new,
            @lpp1_12 = Dvector.new,
            @lCNO0_12 = Dvector.new,
            @lCNO1_12 = Dvector.new,
            @cap_12 = Dvector.new,
            @cond_eff_12 = Dvector.new,
            @ms_Age_12 = Dvector.new,
            @conv_M_s1_12 = Dvector.new,
            @conv_M_e1_12 = Dvector.new,
            @conv_M_s2_12 = Dvector.new,
            @conv_M_e2_12 = Dvector.new,
            @conv_M_s3_12 = Dvector.new,
            @conv_M_e3_12 = Dvector.new ]
        Dvector.read('ellipticalAD_alt2.out', @data12)

	@data13 = [ # columns in ellipticalAD_alt3.out
            @orbit_13 = Dvector.new,
	    @e_13 = Dvector.new,
            @mass_13 = Dvector.new,
            @rhowimp_13 = Dvector.new,
            @v_star_13 = Dvector.new,
            @t1_13 = Dvector.new,
            @tcmin_13 = Dvector.new,
            @lwimp0_13 = Dvector.new,
            @lwimp1_13 = Dvector.new,
            @lpp0_13 = Dvector.new,
            @lpp1_13 = Dvector.new,
            @lCNO0_13 = Dvector.new,
            @lCNO1_13 = Dvector.new,
            @cap_13 = Dvector.new,
            @cond_eff_13 = Dvector.new,
            @ms_Age_13 = Dvector.new,
            @conv_M_s1_13 = Dvector.new,
            @conv_M_e1_13 = Dvector.new,
            @conv_M_s2_13 = Dvector.new,
            @conv_M_e2_13 = Dvector.new,
            @conv_M_s3_13 = Dvector.new,
            @conv_M_e3_13 = Dvector.new ]
        Dvector.read('ellipticalAD_alt3.out', @data13)

	@data14 = [ # columns in etrans.out
            @radius14 = Dvector.new,
	    @etrans14 = Dvector.new ]
        Dvector.read('etrans.out', @data14)

	@etrans14 = @etrans14 * 4.1391085547736407 * 10**(-39)
        @margin = 0.01       
        t.def_enter_page_function { enter_page }
        read_ZAMS('../star_history/star_data')
        
        @colours = {}
        #@colours[0.3] = SeaGreen
        @colours[0.4] = LightGreen
        #@colours[0.5] = DarkBlue
        @colours[0.6] = BrightBlue
        #@colours[0.7] = Periwinkle
        @colours[0.8] = Purple
        #@colours[0.9] = LightPlum
        @colours[1.0] = Burgundy
        @colours[1.2] = DeepPink        
        @colours[1.4] = Red        
        @colours[1.5] = Orange        
	@colours[1.6] = DarkOrange
        @colours[1.8] = Goldenrod
        @colours[2.0] = DarkChocolate

	@linetypes_2 = {}
	@linetypes_2[7] = Line_Type_Solid
	@linetypes_2[8.5] = Line_Type_Dash
	@linetypes_2[10] = Line_Type_Dot

	@linetypes_3 = {}
	@linetypes_3[8] = Line_Type_Solid
	@linetypes_3[9] = Line_Type_Dash
	@linetypes_3[10] = Line_Type_Dot

    end

    def enter_page
        t.page_setup(11*72/2,8.5*72/2)
        t.set_frame_sides(0.15,0.85,0.85,0.15) # left, right, top, bottom in page coords        
    end

    def plot_etrans

        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('WIMP energy transport $\epsilon_\mathrm{trans}$ (erg\,g$^{-1}$\,s$^{-1}$)')
        t.show_xlabel('Height in star, $\log_{10}\big(\frac{r}{R_\star}\big)$')
        
        ydata = etrans14
        xdata = radius14.safe_log10

        xmax = xdata.max; xmin = xdata.min; xmargin = @margin*(xmax-xmin)
        right = -1.0; left = xmin - 3*xmargin 
        ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        top = ymax + 10*ymargin; bottom = ymin
        t.show_plot('left_boundary' => left, 'right_boundary' => right,
         'top_boundary' => top, 'bottom_boundary' => bottom) do
            background
            add_info_line([-2, 0],[0,0])
            t.line_type = Line_Type_Solid
            t.show_polyline(xdata, ydata, Black)
        end       
    end

    def plot_cap

        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('Initial capture rate, $\log_{10}\big(\frac{C(0)}{n\,\mathrm{yr}^{-1}}\big)$')
        t.show_text('text' => '\textbf{Z ='+Metallicity.to_s+'}', 'color' => Black, 'x' => 0.85, 'y' => 0.1)

        xdata1 = (rhowimp*1e-38/270).safe_log10 - 2
	xdata2 = rhowimp.safe_log10
	ydata = cap.safe_log10
	    t.subplot do
		background	        
		t.xaxis_loc = t.xlabel_side = TOP;
		t.left_edge_type = AXIS_HIDDEN        	
		t.bottom_edge_type = AXIS_HIDDEN
		t.right_edge_type = AXIS_HIDDEN
		t.yaxis_type = AXIS_HIDDEN		
		xmax = xdata2.max; xmin = xdata2.min; xmargin = @margin*(xmax-xmin)
        	right = xmax + xmargin; left = xmin - xmargin 
        	ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        	top = ymax + ymargin; bottom = ymin - ymargin
		top = 45.2
		t.show_xlabel('$\log_{10}\big(\frac{\rho_\chi}{\mathrm{GeV}\,\mathrm{cm}^{-3}}\big) = 1.56+ \log_{10}\big(\frac{\rho_\chi}{\mathrm{M}_\odot\,\mathrm{pc}^{-3}}\big)$')
        	t.show_plot('left_boundary' => left, 'right_boundary' => right, 'top_boundary' => top, 'bottom_boundary' => bottom) do end
	    end

        xmax = xdata1.max; xmin = xdata1.min; xmargin = @margin*(xmax-xmin)
        right = xmax + xmargin; left = xmin - xmargin 
        ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
	top = 45.2
        t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.05,'plot_right_margin' => 0,
            'legend_top_margin' => 0.02) do
	    t.subplot do
		t.top_edge_type = AXIS_HIDDEN	
		t.show_plot('left_boundary' => left, 'right_boundary' => right, 'top_boundary' => top, 'bottom_boundary' => bottom) do
		        t.show_xlabel('$\log_{10}\big(\frac{\rho\sigma_\mathrm{SD}}{\bar v}\cdot\mathrm{GeV}^{-1}\,\mathrm{s}^{-1}\,\mathrm{cm}^{2}\big)$')
			colours.sort.reverse.each do |m, col| 
                		sorted_data = [Array.new, Array.new, Array.new]
                		sorted_data[0] = (0..mass.length-1).find_all{|i| (mass[i]/m-1).abs < Epsilon && metal[i] == Metallicity}
                		sorted_data[0].each_index do |i|
                			sorted_data[1][i] = xdata1[sorted_data[0][i]]
                			sorted_data[2][i] = ydata[sorted_data[0][i]]
                			#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
				end
                    		t.show_polyline(sorted_data[1], sorted_data[2], col, m.to_s+'\,M$_\odot$') if sorted_data[1][0] != nil
                	end
		end
	    end
        end       

    end

    def plot_Lwimp
        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')
        t.show_xlabel('$\log_{10}\big(\frac{\rho\sigma_\mathrm{SD}}{\bar v}\cdot\mathrm{GeV}^{-1}\,\mathrm{s}^{-1}\,\mathrm{cm}^{2}\big)$')
        t.show_text('text' => '\textbf{Z ='+Metallicity.to_s+'}', 'color' => Black, 'x' => 0.85, 'y' => 0.1)

        xdata = (rhowimp*1e-38/270).safe_log10 - 2
	ydata = (lwimp0/(lpp0+lCNO0)).safe_log10
	xdata2 = rhowimp.safe_log10
	t.subplot do
		background	        
		t.xaxis_loc = t.xlabel_side = TOP;
		t.left_edge_type = AXIS_HIDDEN        	
		t.bottom_edge_type = AXIS_HIDDEN
		t.right_edge_type = AXIS_HIDDEN		
		t.yaxis_type = AXIS_HIDDEN
		xmax = xdata2.max; xmin = xdata2.min; xmargin = @margin*(xmax-xmin)
        	right = xmax + xmargin; left = xmin - xmargin 
        	ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        	top = ymax + ymargin; bottom = ymin - ymargin
		top = 4		
		t.show_xlabel('$\log_{10}\big(\frac{\rho_\chi}{\mathrm{GeV}\,\mathrm{cm}^{-3}}\big) = 1.56+ \log_{10}\big(\frac{\rho_\chi}{\mathrm{M}_\odot\,\mathrm{pc}^{-3}}\big)$')
        	t.show_plot('left_boundary' => left, 'right_boundary' => right, 'top_boundary' => top, 'bottom_boundary' => bottom) do end
	end

        xmax = xdata.max; xmin = xdata.min; xmargin = @margin*(xmax-xmin)
        right = xmax + xmargin; left = xmin - xmargin 
        ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
	top = 4        
	t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.05,'plot_right_margin' => 0,
	 'legend_top_margin' => 0.02) do
            t.top_edge_type = AXIS_HIDDEN	
	    t.show_plot('left_boundary' => left, 'right_boundary' => right,
             'top_boundary' => top, 'bottom_boundary' => bottom) do
                colours.sort.reverse.each do |m, col| 
                    sorted_data = [Array.new, Array.new, Array.new]
                    sorted_data[0] = (0..mass.length-1).find_all{|i| (mass[i]/m-1).abs < Epsilon && metal[i] == Metallicity}
                    sorted_data[0].each_index do |i|
                        sorted_data[1][i] = xdata[sorted_data[0][i]]
                        sorted_data[2][i] = ydata[sorted_data[0][i]]
                	#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
			#print m, "   ", sorted_data[1][i]+40+270.safe_log10,"   ",sorted_data[2][i],"\n"
                    end
                    t.show_polyline(sorted_data[1], sorted_data[2], col, m.to_s+'\,M$_\odot$') if sorted_data[1][0] != nil
                end
            end
        end       
    end


    def plot_ann
        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')
        t.show_xlabel('Annihilation cross-section, $\log_{10}\big(\frac{\langle\sigma_\mathrm{a} v\rangle_0}{\mathrm{cm}^3\,\mathrm{s}^{-1}}\big)$')
        t.show_text('text' => '\textbf{Z = 0.01}', 'color' => Black, 'x' => 0.2, 'y' => 0.77)
        t.show_text('text' => '\textbf{M = 1.0 $\mathbf{M_{\odot}}$}', 'color' => Black, 'x' => 0.2, 'y' => 0.87)

        ydata = (lwimp0_3/(lpp0_3+lCNO0_3)).safe_log10
        xdata = sigmav_3.safe_log10

        xmax = xdata.max; xmin = xdata.min; xmargin = @margin*(xmax-xmin)
        right = xmax + xmargin; left = 3e-27.safe_log10 
        ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        top = 3; bottom = -1.8
        t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.45,'plot_right_margin' => 0,
         'legend_top_margin' => 0.02) do
            t.show_plot('left_boundary' => left, 'right_boundary' => right,
             'top_boundary' => top, 'bottom_boundary' => bottom) do
                background
                linetypes_3.sort.reverse.each do |rho, linetype|
		    t.line_type = linetype
                    sorted_data = [Array.new, Array.new, Array.new]
                    sorted_data[0] = (0..mass_3.length-1).find_all{|i| (rhowimp_3[i].safe_log10/rho-1).abs < Epsilon}
                    sorted_data[0].each_index do |i|
                        sorted_data[1][i] = xdata[sorted_data[0][i]]
                        sorted_data[2][i] = ydata[sorted_data[0][i]]
                	#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
			#print sorted_data[1][i]+40+v_star[0].safe_log10,"   ",sorted_data[2][i],"\n"
                    end
		    if sorted_data[1][0] != nil then
			x = sorted_data[1].sort
			y = sorted_data[2].sort {|a,b| sorted_data[1][sorted_data[2].index(a)] <=> sorted_data[1][sorted_data[2].index(b)]}
			t.show_polyline(x, y, Black, '$\log_{10}\big(\rho_\chi\cdot\mathrm{GeV}^{-1}\,\mathrm{cm}^{3}\big) = '+rho.to_s+'$')
		    end
                end
            t.line_type = Line_Type_Solid
	    end
        end       
    end


    def plot_Mwimp
        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')
        t.show_xlabel('WIMP mass, $\log_{10}\big(\frac{m_\chi}{\mathrm{GeV}}\big)$')
        t.show_text('text' => '\textbf{Z = 0.01}', 'color' => Black, 'x' => 0.2, 'y' => 0.1)
        t.show_text('text' => '\textbf{M = 1.0 $\mathbf{M_{\odot}}$}', 'color' => Black, 'x' => 0.2, 'y' => 0.2)

        ydata = (lwimp0_2/(lpp0_2+lCNO0_2)).safe_log10
        xdata = mchi_2.safe_log10

        xmax = xdata.max; xmin = xdata.min; xmargin = @margin*(xmax-xmin)
        right = xmax + xmargin; left = 40.safe_log10 
        ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
        t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.45,'plot_right_margin' => 0,
         'legend_top_margin' => 0.02) do
            t.show_plot('left_boundary' => left, 'right_boundary' => right,
             'top_boundary' => top, 'bottom_boundary' => bottom) do
                background
                linetypes_2.sort.reverse.each do |rho, linetype|
		    t.line_type = linetype
                    sorted_data = [Array.new, Array.new, Array.new]
                    sorted_data[0] = (0..mass_2.length-1).find_all{|i| (rhowimp_2[i].safe_log10/rho-1).abs < Epsilon}
                    sorted_data[0].each_index do |i|
                        sorted_data[1][i] = xdata[sorted_data[0][i]]
                        sorted_data[2][i] = ydata[sorted_data[0][i]]
                	#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
			#print sorted_data[1][i]+40+v_star[0].safe_log10,"   ",sorted_data[2][i],"\n"
                    end
		    if sorted_data[1][0] != nil then
			x = sorted_data[1].sort
			y = sorted_data[2].sort {|a,b| sorted_data[1][sorted_data[2].index(a)] <=> sorted_data[1][sorted_data[2].index(b)]}
			t.show_polyline(x, y, Black, '$\log_{10}\big(\rho_\chi\cdot\mathrm{GeV}^{-1}\,\mathrm{cm}^{3}\big) = '+rho.to_s+'$')
		    end
                end
            t.line_type = Line_Type_Solid
	    end
        end       
    end


    def plot_HRPlotp1
	stars = Array.new
	relevant_indexes = (0..mass_4.length-1).find_all{|i| ((lwimp0_4[i]/(lpp0_4[i]+lCNO0_4[i])).safe_log10 - 1).abs < HREpsilon}
	relevant_indexes.each_index do |i| stars[i] = hrdata[relevant_indexes[i]] end
	background
	t.show_text('text' => '0\,yr', 'color' => Black, 'x' => 0.27, 'y' => 0.735, 'scale'=>0.7)
	t.show_text('text' => '0.5\,Myr', 'color' => Black, 'x' => 0.35, 'y' => 0.64, 'scale'=>0.7)
	t.show_text('text' => '1.0\,Myr', 'color' => Black, 'x' => 0.52, 'y' => 0.585, 'scale'=>0.7)
	t.show_text('text' => '2.0\,Myr', 'color' => Black, 'x' => 0.595, 'y' => 0.463, 'scale'=>0.7)
	t.show_text('text' => '5.0\,Myr', 'color' => Black, 'x' => 0.708, 'y' => 0.508, 'scale'=>0.7)
	t.show_text('text' => '$\ge 44$\,Myr', 'color' => Black, 'x' => 0.675, 'y' => 0.623, 'scale'=>0.7)
	individ_HR(stars,1)
    end

    def plot_Trhop1
	stars = Array.new
	relevant_indexes = (0..mass_4.length-1).find_all{|i| ((lwimp0_4[i]/(lpp0_4[i]+lCNO0_4[i])).safe_log10 - 1).abs < HREpsilon}
	relevant_indexes.each_index do |i| stars[i] = hrdata[relevant_indexes[i]] end
	background
	t.show_text('text' => '0\,yr', 'color' => Black, 'x' => 0.95, 'y' => 0.82, 'scale'=>0.7)
	t.show_text('text' => '0.5\,Myr', 'color' => Black, 'x' => 0.85, 'y' => 0.73, 'scale'=>0.7)
	t.show_text('text' => '1.0\,Myr', 'color' => Black, 'x' => 0.77, 'y' => 0.65, 'scale'=>0.7)
	t.show_text('text' => '2.0\,Myr', 'color' => Black, 'x' => 0.64, 'y' => 0.51, 'scale'=>0.7)
	t.show_text('text' => '5.0\,Myr', 'color' => Black, 'x' => 0.39, 'y' => 0.25, 'scale'=>0.7)
	t.show_text('text' => '$\ge 44$\,Myr', 'color' => Black, 'x' => 0.15, 'y' => 0.06, 'scale'=>0.7)
	individ_Trho(stars, 1, textloc=4, legloc=1)
    end

    def plot_HRPlot0
	stars = Array.new
	relevant_indexes = (0..mass_4.length-1).find_all{|i| ((lwimp0_4[i]/(lpp0_4[i]+lCNO0_4[i])).safe_log10).abs < HREpsilon}
	relevant_indexes.each_index do |i| stars[i] = hrdata[relevant_indexes[i]] end
	background
	t.show_text('text' => '0\,yr', 'color' => Black, 'x' => 0.285, 'y' => 0.62, 'scale'=>0.7)
	t.show_text('text' => '1.0\,Myr', 'color' => Black, 'x' => 0.35, 'y' => 0.5, 'scale'=>0.7)
	t.show_text('text' => '3.3\,Myr', 'color' => Black, 'x' => 0.44, 'y' => 0.62, 'scale'=>0.72)
	t.show_text('text' => '2.0\,Gyr', 'color' => Black, 'x' => 0.3, 'y' => 0.75, 'scale'=>0.7)
	t.show_text('text' => '3.4\,Gyr', 'color' => Black, 'x' => 0.457, 'y' => 0.668, 'scale'=>0.7)
	t.show_text('text' => '4.05\,Gyr', 'color' => Black, 'x' => 0.55, 'y' => 0.7, 'scale'=>0.7)
	t.show_text('text' => '4.11\,Gyr', 'color' => Black, 'x' => 0.41, 'y' => 0.778, 'scale'=>0.7)
	t.show_text('text' => '4.30\,Gyr', 'color' => Black, 'x' => 0.625, 'y' => 0.745, 'scale'=>0.7)
	t.show_arrow('head' => [0.338,0.678], 'tail'=> [0.3, 0.74], 'head_scale' => 0.5, 'line_width' => 0.6, 'tail_marker' => 'None')
	t.show_arrow('head' => [0.357,0.6], 'tail'=> [0.35,0.53], 'head_scale' => 0.5, 'line_width' => 0.6, 'tail_marker' => 'None')
	individ_HR(stars, 0)
    end

    def plot_Trho0
	stars = Array.new
	relevant_indexes = (0..mass_4.length-1).find_all{|i| ((lwimp0_4[i]/(lpp0_4[i]+lCNO0_4[i])).safe_log10).abs < HREpsilon}
	relevant_indexes.each_index do |i| stars[i] = hrdata[relevant_indexes[i]] end
	background
	t.show_text('text' => '0\,yr', 'color' => Black, 'x' => 0.205, 'y' => 0.425, 'scale'=>0.7)
	t.show_text('text' => '1.0\,Myr', 'color' => Black, 'x' => 0.07, 'y' => 0.475, 'scale'=>0.7)
	t.show_text('text' => '3.3\,Myr', 'color' => Black, 'x' => 0.1, 'y' => 0.375, 'scale'=>0.7)
	t.show_text('text' => '2.0\,Gyr', 'color' => Black, 'x' => 0.24, 'y' => 0.47, 'scale'=>0.7)
	t.show_text('text' => '3.4\,Gyr', 'color' => Black, 'x' => 0.29, 'y' => 0.53, 'scale'=>0.7)
	t.show_text('text' => '4.05\,Gyr', 'color' => Black, 'x' => 0.38, 'y' => 0.62, 'scale'=>0.7)
	t.show_text('text' => '4.11\,Gyr', 'color' => Black, 'x' => 0.5, 'y' => 0.7, 'scale'=>0.7)
	t.show_text('text' => '4.30\,Gyr', 'color' => Black, 'x' => 0.92, 'y' => 0.86, 'scale'=>0.7)
	t.show_arrow('head' => [0.123,0.447], 'tail'=> [0.075, 0.465], 'head_scale' => 0.5, 'line_width' => 0.6, 'tail_marker' => 'None')
	individ_Trho(stars, 0, textloc=1, legloc=4)
    end

    def plot_HRPlotn1
	stars = Array.new
	relevant_indexes = (0..mass_4.length-1).find_all{|i| ((lwimp0_4[i]/(lpp0_4[i]+lCNO0_4[i])).safe_log10 + 1).abs < HREpsilon}
	relevant_indexes.each_index do |i| stars[i] = hrdata[relevant_indexes[i]] end
	background
	t.show_text('text' => '0\,yr', 'color' => Black, 'x' => 0.35, 'y' => 0.615, 'scale'=>0.7)
	t.show_text('text' => '1.0\,Gyr', 'color' => Black, 'x' => 0.2, 'y' => 0.66, 'scale'=>0.7)
	t.show_text('text' => '2.1\,Gyr', 'color' => Black, 'x' => 0.42, 'y' => 0.665, 'scale'=>0.7)
	t.show_text('text' => '2.9\,Gyr', 'color' => Black, 'x' => 0.54, 'y' => 0.7, 'scale'=>0.7)
	t.show_text('text' => '3.0\,Gyr', 'color' => Black, 'x' => 0.45, 'y' => 0.787, 'scale'=>0.7)
	t.show_text('text' => '3.1\,Gyr', 'color' => Black, 'x' => 0.625, 'y' => 0.745, 'scale'=>0.7)
	t.show_arrow('head' => [0.31,0.675], 'tail'=> [0.245, 0.67], 'head_scale' => 0.5, 'line_width' => 0.6, 'tail_marker' => 'None')
	individ_HR(stars, -1)
    end

    def plot_Trhon1
	stars = Array.new
	relevant_indexes = (0..mass_4.length-1).find_all{|i| ((lwimp0_4[i]/(lpp0_4[i]+lCNO0_4[i])).safe_log10 + 1).abs < HREpsilon}
	relevant_indexes.each_index do |i| stars[i] = hrdata[relevant_indexes[i]] end
	background
	t.show_text('text' => '0\,yr', 'color' => Black, 'x' => 0.14, 'y' => 0.48, 'scale'=>0.7)
	t.show_text('text' => '1.0\,Gyr', 'color' => Black, 'x' => 0.195, 'y' => 0.54, 'scale'=>0.7)
	t.show_text('text' => '2.1\,Gyr', 'color' => Black, 'x' => 0.225, 'y' => 0.6, 'scale'=>0.7)
	t.show_text('text' => '2.9\,Gyr', 'color' => Black, 'x' => 0.295, 'y' => 0.70, 'scale'=>0.7)
	t.show_text('text' => '3.0\,Gyr', 'color' => Black, 'x' => 0.6, 'y' => 0.85, 'scale'=>0.7)
	t.show_text('text' => '3.1\,Gyr', 'color' => Black, 'x' => 0.835, 'y' => 0.86, 'scale'=>0.7)
	individ_Trho(stars, -1, textloc=4, legloc=4)
    end

    def plot_HRPlotn3
	stars = Array.new
	relevant_indexes = (0..mass_4.length-1).find_all{|i| ((lwimp0_4[i]/(lpp0_4[i]+lCNO0_4[i])).safe_log10 + 3).abs < HREpsilon}
	relevant_indexes.each_index do |i| stars[i] = hrdata[relevant_indexes[i]] end
	background
	t.show_text('text' => '0\,yr', 'color' => Black, 'x' => 0.35, 'y' => 0.615, 'scale'=>0.7)
	t.show_text('text' => '1.0\,Gyr', 'color' => Black, 'x' => 0.2, 'y' => 0.66, 'scale'=>0.7)
	t.show_text('text' => '2.0\,Gyr', 'color' => Black, 'x' => 0.43, 'y' => 0.665, 'scale'=>0.7)
	t.show_text('text' => '2.69\,Gyr', 'color' => Black, 'x' => 0.55, 'y' => 0.7, 'scale'=>0.7)
	t.show_text('text' => '2.77\,Gyr', 'color' => Black, 'x' => 0.45, 'y' => 0.787, 'scale'=>0.7)
	t.show_text('text' => '2.82\,Gyr', 'color' => Black, 'x' => 0.635, 'y' => 0.745, 'scale'=>0.7)
	t.show_arrow('head' => [0.307,0.675], 'tail'=> [0.245, 0.67], 'head_scale' => 0.5, 'line_width' => 0.6, 'tail_marker' => 'None')
	individ_HR(stars, -3)
    end

    def plot_Trhon3
	stars = Array.new
	relevant_indexes = (0..mass_4.length-1).find_all{|i| ((lwimp0_4[i]/(lpp0_4[i]+lCNO0_4[i])).safe_log10 + 3).abs < HREpsilon}
	relevant_indexes.each_index do |i| stars[i] = hrdata[relevant_indexes[i]] end
	background
	t.show_text('text' => '0\,yr', 'color' => Black, 'x' => 0.14, 'y' => 0.49, 'scale'=>0.7)
	t.show_text('text' => '1.0\,Gyr', 'color' => Black, 'x' => 0.19, 'y' => 0.55, 'scale'=>0.7)
	t.show_text('text' => '2.0\,Gyr', 'color' => Black, 'x' => 0.22, 'y' => 0.62, 'scale'=>0.7)
	t.show_text('text' => '2.69\,Gyr', 'color' => Black, 'x' => 0.295, 'y' => 0.71, 'scale'=>0.7)
	t.show_text('text' => '2.77\,Gyr', 'color' => Black, 'x' => 0.42, 'y' => 0.83, 'scale'=>0.7)
	t.show_text('text' => '2.82\,Gyr', 'color' => Black, 'x' => 0.9, 'y' => 0.785, 'scale'=>0.7)
	individ_Trho(stars, -3, textloc=5, legloc=4)
    end


    def individ_HR(d, lw, right_margin=0, left_margin=0)
        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        margin = 0.1
	t.show_xlabel('Effective (surface) temperature, $\log_{10}\big(\frac{T_\mathrm{eff}}{\mathrm{K}}\big)$')
	t.show_ylabel('Luminosity, $\log_{10}(L/\mathrm{L}_\odot)$')
	tmax = 0; lmax = 0; tmin = 1e20; lmin = 1e20
	d.each do |star|
		tmax = [tmax, star[1].max].max
		lmax = [lmax, star[0].max].max
		tmin = [tmin, star[1].min].min
		lmin = [lmin, star[0].min].min
	end
	xmargin = margin*(tmax-tmin); left = tmax + xmargin; right = tmin - xmargin 
	ymargin = margin*(lmax-lmin); top = lmax + ymargin; bottom = lmin - ymargin
	t.show_text('text' => '\textbf{Z = 0.01}', 'color' => Black, 'x' => 0.19, 'y' => 0.1)
	t.show_text('text' => '\textbf{log}$\mathbf{_{10}\Big(\frac{L_\mathbf{W,max}}{L_\mathbf{nuc}(0)}\Big)} \mathbf{\approx} \mathbf{'+lw.to_s+'}$', 'color' => Black, 'x' => 0.19, 'y' => 0.2)

	t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.73,'plot_right_margin' => 0,
         'legend_top_margin' => 0.02) do
		t.show_plot('left_boundary' => left, 'right_boundary' => right,
		 'top_boundary' => top, 'bottom_boundary' => bottom) do
			d.reverse.each do |star|        
			    xs = star[1]
		            ys = star[0]
			    stroke_track(xs, ys, colours[star[4]], star[4].to_s+'\,M$_\odot$', Line_Type_Solid, 0)
			    star[5].each do |i| mark_spot(xs[i-1], ys[i-1], colours[star[4]]) end
			    #dots_in_gap(xs, ys)
		            #mark_points(xs, ys)
	        	end
			add_info_line(@zams_log_Surface_Temp, @zams_log_Luminosity)
			t.save_legend_info('text' => 'ZAMS')			    
	        end
	end
    end
    

    def individ_Trho(d, lw, textloc=3, legloc=2, right_margin=0, left_margin=0)
        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        margin = 0.1
        t.show_xlabel('Central density, $\log_{10}\big(\frac{\rho_\mathrm{c}}{\mathrm{g}\,\mathrm{cm}^{-3}}\big)$')
        t.show_ylabel('Central temperature, $\log_{10}\big(\frac{T_\mathrm{c}}{\mathrm{K}}\big)$')
	
	case textloc	
		when 1 then
 			textx = 0.19; texty = 0.9; texty2 = 0.82
		when 2 then
 			textx = 0.81; texty = 0.9; texty2 = 0.82
		when 3 then
 			textx = 0.19; texty = 0.1; texty2 = 0.2
		when 4 then
 			textx = 0.81; texty = 0.1; texty2 = 0.2
		when 5 then
			textx = 0.50; texty = 0.1; texty2 = 0.2
	end

	case legloc	
		when 1 then
 			legx = 0.05; legy = 0.02
		when 2 then
 			legx = 0.73; legy = 0.02
		when 3 then
 			legx = 0.05; legy = 0.55
		when 4 then
			if (textloc == 4) then
				legx = 0.71; legy = 0.28
			else
	 			legx = 0.73; legy = 0.52
			end
	end

	rmax = 0; tmax = 0; rmin = 1e20; tmin = 1e20
	d.each do |star|
		rmax = [rmax, star[3].max].max
		tmax = [tmax, star[2].max].max
		rmin = [rmin, star[3].min].min
		tmin = [tmin, star[2].min].min
	end
	xmargin = margin*(rmax-rmin); right = rmax + xmargin; left = rmin - xmargin 
	ymargin = margin*(tmax-tmin); top = tmax + ymargin; bottom = tmin - ymargin
	t.show_text('text' => '\textbf{Z = 0.01}', 'color' => Black, 'x' => textx, 'y' => texty)
	t.show_text('text' => '\textbf{log}$\mathbf{_{10}\Big(\frac{L_\mathbf{W,max}}{L_\mathbf{nuc}(0)}\Big)} \mathbf{\approx} \mathbf{'+lw.to_s+'}$', 'color' => Black, 'x' => textx, 'y' => texty2)

	t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => legx,'plot_right_margin' => 0,
         'legend_top_margin' => legy) do
		t.show_plot('left_boundary' => left, 'right_boundary' => right,
		 'top_boundary' => top, 'bottom_boundary' => bottom) do
			d.reverse.each do |star|        
			    xs = star[3]
		            ys = star[2]
		            stroke_track(xs, ys, colours[star[4]], star[4].to_s+'\,M$_\odot$', Line_Type_Solid, 0)
		            star[5].each do |i| mark_spot(xs[i-1], ys[i-1], colours[star[4]]) end
			    #dots_in_gap(xs, ys)
		            #mark_points(xs, ys)
			end
			add_info_line(@zams_log_Center_Density, @zams_log_Center_Temp)
			t.save_legend_info('text' => 'ZAMS')    
	        end
	end
    end


    def plot_pp
        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('pp--chain luminosity, $\log_{10}\Big[\frac{L_\mathrm{pp}(t_\mathrm{adjust})}{L_\mathrm{pp}(0)}\Big]$')
        t.show_xlabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')
        t.show_text('text' => '\textbf{Z ='+Metallicity.to_s+'}', 'color' => Black, 'x' => 0.15, 'y' => 0.1)

        xdata = (lwimp0/(lpp0+lCNO0)).safe_log10
        ydata = (lpp1/lpp0).safe_log10

        xmax = xdata.max; xmin = xdata.min; xmargin = @margin*(xmax-xmin)
        right = xmax + xmargin; left = xmin - xmargin 
        ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
        t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.73,'plot_right_margin' => 0,
            'legend_top_margin' => 0.02) do
            t.show_plot('left_boundary' => left, 'right_boundary' => right,
                'top_boundary' => top, 'bottom_boundary' => bottom) do
                background
                colours.sort.reverse.each do |m, col| 
                    sorted_data = [Array.new, Array.new, Array.new]
                    sorted_data[0] = (0..mass.length-1).find_all{|i| (mass[i]/m-1).abs < Epsilon && metal[i] == Metallicity}
                    sorted_data[0].each_index do |i|
                        sorted_data[1][i] = xdata[sorted_data[0][i]]
                        sorted_data[2][i] = ydata[sorted_data[0][i]]
                	#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
			#print sorted_data[1][i]+40+v_star[0].safe_log10,"   ",sorted_data[2][i],"\n"
                    end
                    t.show_polyline(sorted_data[1], sorted_data[2], col, m.to_s+'\,M$_\odot$') if sorted_data[1][0] != nil
                end
            end
        end       
    end


    def plot_CNO
        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('CNO--process luminosity, $\log_{10}\Big[\frac{L_\mathrm{CNO}(t_\mathrm{adjust})}{L_\mathrm{CNO}(0)}\Big]$')
        t.show_xlabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')
        t.show_text('text' => '\textbf{Z ='+Metallicity.to_s+'}', 'color' => Black, 'x' => 0.15, 'y' => 0.1)

        xdata = (lwimp0/(lpp0+lCNO0)).safe_log10
        ydata = (lCNO1/lCNO0).safe_log10

        xmax = xdata.max; xmin = xdata.min; xmargin = @margin*(xmax-xmin)
        right = xmax + xmargin; left = xmin - xmargin 
        ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
        t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.73,'plot_right_margin' => 0,
            'legend_top_margin' => 0.02) do
            t.show_plot('left_boundary' => left, 'right_boundary' => right,
                'top_boundary' => top, 'bottom_boundary' => bottom) do
                background
                colours.sort.reverse.each do |m, col| 
                    sorted_data = [Array.new, Array.new, Array.new]
                    sorted_data[0] = (0..mass.length-1).find_all{|i| (mass[i]/m-1).abs < Epsilon && metal[i] == Metallicity}
                    sorted_data[0].each_index do |i|
                        sorted_data[1][i] = xdata[sorted_data[0][i]]
                        sorted_data[2][i] = ydata[sorted_data[0][i]]
                	#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
			#print sorted_data[1][i]+40+v_star[0].safe_log10,"   ",sorted_data[2][i],"\n"
                    end
                    t.show_polyline(sorted_data[1], sorted_data[2], col, m.to_s+'\,M$_\odot$') if sorted_data[1][0] != nil
                end
            end
        end       
    end


    def plot_Tcmin
        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('Central temperature, $\log_{10}\Big[\frac{T_\mathrm{c}(t_\mathrm{adjust})}{\mathrm{K}}\Big]$')
        t.show_xlabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')
        t.show_text('text' => '\textbf{Z ='+Metallicity.to_s+'}', 'color' => Black, 'x' => 0.15, 'y' => 0.1)

        xdata = (lwimp0/(lpp0+lCNO0)).safe_log10
        ydata = tcmin

        xmax = xdata.max; xmin = xdata.min; xmargin = @margin*(xmax-xmin)
        right = xmax + xmargin; left = xmin - xmargin 
        ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
        t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.73,'plot_right_margin' => 0,
            'legend_top_margin' => 0.02) do
            t.show_plot('left_boundary' => left, 'right_boundary' => right,
                'top_boundary' => top, 'bottom_boundary' => bottom) do
                background
                colours.sort.reverse.each do |m, col| 
                    sorted_data = [Array.new, Array.new, Array.new]
                    sorted_data[0] = (0..mass.length-1).find_all{|i| (mass[i]/m-1).abs < Epsilon && metal[i] == Metallicity}
                    sorted_data[0].each_index do |i|
                        sorted_data[1][i] = xdata[sorted_data[0][i]]
                        sorted_data[2][i] = ydata[sorted_data[0][i]]
                	#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
			#print m, "   ", sorted_data[1][i],"   ",sorted_data[2][i],"\n"
                    end
                    t.show_polyline(sorted_data[1], sorted_data[2], col, m.to_s+'\,M$_\odot$') if sorted_data[1][0] != nil
                end
            end
        end       
    end


    def plot_cond
        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('WIMP conductive effectiveness, $\log_{10}\big[\mathfrak{E}(t_\mathrm{adjust})\big]$')
        t.show_xlabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')
        t.show_text('text' => '\textbf{Z ='+Metallicity.to_s+'}', 'color' => Black, 'x' => 0.85, 'y' => 0.1)

        xdata = (lwimp0/(lpp0+lCNO0)).safe_log10
        ydata = cond_eff.safe_log10 - 5.3830931833791515

        xmax = 2; xmin = -5; xmargin = @margin*(xmax-xmin)
        right = xmax + xmargin; left = xmin - xmargin 
        ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
        t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.05,'plot_right_margin' => 0,
            'legend_top_margin' => 0.02) do
            t.show_plot('left_boundary' => left, 'right_boundary' => right,
                'top_boundary' => top, 'bottom_boundary' => bottom) do
                background
                colours.sort.reverse.each do |m, col| 
                    sorted_data = [Array.new, Array.new, Array.new]
                    sorted_data[0] = (0..mass.length-1).find_all{|i| (mass[i]/m-1).abs < Epsilon && metal[i] == Metallicity}
                    sorted_data[0].each_index do |i|
                        sorted_data[1][i] = xdata[sorted_data[0][i]]
                        sorted_data[2][i] = ydata[sorted_data[0][i]]
                	#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
			#print m, "   ", sorted_data[1][i],"   ",sorted_data[2][i],"\n"
                    end
                    t.show_polyline(sorted_data[1], sorted_data[2], col, m.to_s+'\,M$_\odot$') if sorted_data[1][0] != nil
                end
            end
        end       
    end


    def plot_conv03
	individ_conv(0.3)
    end

    def plot_conv04
	individ_conv(0.4)
    end

    def plot_conv05
	individ_conv(0.5)
    end

    def plot_conv06
	individ_conv(0.6)
    end

    def plot_conv07
	individ_conv(0.7)
    end

    def plot_conv08
	individ_conv(0.8)
    end

    def plot_conv09
	individ_conv(0.9)
    end

    def plot_conv10
	individ_conv(1.0)
    end

    def plot_conv12
	individ_conv(1.2)
    end

    def plot_conv14
	individ_conv(1.4)
    end

    def plot_conv16
	individ_conv(1.6)
    end

    def plot_conv18
	individ_conv(1.8)
    end

    def plot_conv20
	individ_conv(2.0)
    end

    def individ_conv(m)

        t.ylabel_scale = 2
        t.xlabel_scale = 2
	t.set_frame_sides(0.15,0.95,0.95,0.2)
        t.yaxis_numeric_label_scale = 1.8
        t.xaxis_numeric_label_scale = 1.8
        t.show_ylabel('Enclosed mass ($M/\mathrm{M}_\odot$)')
        t.show_xlabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')	
	if m > 1.4 
		offset = 0.1 
	else
		offset = 0.0
	end
        t.show_text('text' => '\textbf{M = '+m.to_s+'\,$\mathbf{M_{\odot}}$}', 'color' => Black, 'x' => 0.3, 'y' => 0.23+offset, 'scale' => 1.8)
	t.show_text('text' => '\textbf{Z = '+Metallicity.to_s+'}', 'color' => Black, 'x' => 0.3, 'y' => 0.1+offset, 'scale' => 1.8)
	col = colours[m]

        xdata = (lwimp0/(lpp0+lCNO0)).safe_log10
	sorted_data = [Array.new, Array.new, Array.new, Array.new, Array.new, Array.new, Array.new, Array.new]
        sorted_data[0] = (0..mass.length-1).find_all{|i| (mass[i]/m-1).abs < Epsilon && metal[i] == Metallicity}
        sorted_data[0].each_index do |i|
		sorted_data[1][i] = xdata[sorted_data[0][i]]
		sorted_data[2][i] = conv_M_s1[sorted_data[0][i]]
		sorted_data[3][i] = conv_M_e1[sorted_data[0][i]]
		sorted_data[4][i] = conv_M_s2[sorted_data[0][i]]
		sorted_data[5][i] = conv_M_e2[sorted_data[0][i]]
		sorted_data[6][i] = conv_M_s3[sorted_data[0][i]]
		sorted_data[7][i] = conv_M_e3[sorted_data[0][i]]
	end

        ary = sorted_data[2..7]        	
	ymax = max_of_many(ary) * 1.05
        xmax = sorted_data[1].max; xmin = sorted_data[1].min
        
	if m == 0.8 
		t.yaxis_locations_for_major_ticks = [0.2, 0.4, 0.6, 0.8]
	end
	t.xaxis_type = AXIS_WITH_MAJOR_TICKS_AND_NUMERIC_LABELS
	t.yaxis_type = AXIS_WITH_MAJOR_TICKS_AND_NUMERIC_LABELS
	t.top_edge_type = AXIS_WITH_MAJOR_TICKS_ONLY
	t.right_edge_type = AXIS_WITH_MAJOR_TICKS_ONLY

        t.show_plot('boundaries' => [ xmin, xmax, ymax, 0.001 ]) do
            background
            draw_conv_info(sorted_data[1], sorted_data[2], sorted_data[3], sorted_data[4], sorted_data[5], sorted_data[6], sorted_data[7], col, ymax)
        end

    end


    def draw_conv_info(xs, s1, e1, s2, e2, s3, e3, col, ymax)

        t.stroke_color = col
	t.fill_color = col
        t.stroke_width = 2
        null_zone = -20
	zones = [[Array.new, Array.new], [Array.new, Array.new], [Array.new, Array.new]]

	xs.each_index do |i| 

		temp1 = [s1[i], s2[i], s3[i]].find_all{|j| j != null_zone} 
		temp2 = [e1[i], e2[i], e3[i]].find_all{|j| j != null_zone}
		
		case temp1.size	
			when 0 then
 				zones[0][0][i] = null_zone
				zones[0][1][i] = null_zone
				zones[1][0][i] = null_zone
				zones[1][1][i] = null_zone
				zones[2][0][i] = null_zone
				zones[2][1][i] = null_zone
			when 1 then
				if (temp1[0] < Epsilon) then
					zones[0][0][i] = temp1[0]
					zones[0][1][i] = temp2[0]
					zones[1][0][i] = null_zone
					zones[1][1][i] = null_zone
					if (temp2[0]/ymax/1.05 - 1 < Epsilon) then
						zones[2][0][i] = temp1[0]
						zones[2][1][i] = temp2[0]
					else
						zones[2][0][i] = null_zone
						zones[2][1][i] = null_zone
					end
				else
					zones[0][0][i] = null_zone
					zones[0][1][i] = null_zone
					zones[1][0][i] = null_zone
					zones[1][1][i] = null_zone
					zones[2][0][i] = temp1[0]
					zones[2][1][i] = temp2[0]
				end
			when 2 then
				zones[0][0][i] = temp1.min
				zones[0][1][i] = temp2.min
				zones[1][0][i] = null_zone
				zones[1][1][i] = null_zone
				zones[2][0][i] = temp1.max
				zones[2][1][i] = temp2.max
			when 3 then
				zones[0][0][i] = temp1.min
				zones[0][1][i] = temp2.min
				zones[2][0][i] = temp1.max
				zones[2][1][i] = temp2.max
				zones[1][0][i] = temp1.detect{|j| j != zones[0][0][i] && j != zones[2][0][i]}
				zones[1][1][i] = temp2.detect{|j| j != zones[0][1][i] && j != zones[2][1][i]}
		end

	end

	[0,1,2].each_index do |j|
		magicspot = (0..xs.length-1).detect{|i| zones[j][1][i] != null_zone}
		if (magicspot != nil) then 	
			t.move_to_point(xs[magicspot], zones[j][1][magicspot])
			xs.each_index do |i| 
				if (zones[j][1][i] != null_zone) then t.append_point_to_path(xs[i], zones[j][1][i]) end 
			end
			xs.reverse.each_with_index do |x,i| 
				temp = xs.length-1-i
				if (zones[j][0][temp] != null_zone) then t.append_point_to_path(x, zones[j][0][temp]) end
			end
			t.close_fill_and_stroke
		        t.clip
		end
	end

    end


    def plot_TMS
        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('Main-sequence lifetime (Gyr)')
        t.show_xlabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')
        t.show_text('text' => '\textbf{Z ='+Metallicity.to_s+'}', 'color' => Black, 'x' => 0.15, 'y' => 0.88)

        xdata = (lwimp0/(lpp0+lCNO0)).safe_log10
	ydata = ms_Age/1e9
        ydata.each_index do |i| ydata[i] = 20 if ydata[i] < 0 end

        xmax = 1.5; xmin = -6.5; xmargin = @margin*(xmax-xmin)
        right = xmax + xmargin; left = xmin - xmargin 
        ymax = 14; ymin = 0; ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
        t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.05,'plot_right_margin' => 0,
            'legend_top_margin' => 0.13) do
            t.show_plot('left_boundary' => left, 'right_boundary' => right,
                'top_boundary' => top, 'bottom_boundary' => bottom) do
                background
	        t.line_type = Line_Type_Dash
	        t.line_color = SlateGray
	        t.line_width = 1
	        t.stroke_line(-100, 13.7, 100, 13.7)
		t.line_type = Line_Type_Solid       
                colours.sort.reverse.each do |m, col| 
                    sorted_data = [Array.new, Array.new, Array.new]
                    sorted_data[0] = (0..mass.length-1).find_all{|i| (mass[i]/m-1).abs < Epsilon && metal[i] == Metallicity}
                    sorted_data[0].each_index do |i|
                        sorted_data[1][i] = xdata[sorted_data[0][i]]
                        sorted_data[2][i] = ydata[sorted_data[0][i]]
                        #mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
                        #print  m, "   ", sorted_data[1][i],"   ",sorted_data[2][i], "   ",rhowimp[sorted_data[0][i]].safe_log10,"\n"
                    end
                    t.show_polyline(sorted_data[1], sorted_data[2], col, m.to_s+'\,M$_\odot$') if sorted_data[1][0] != nil
                end
            end
        end
    end


    def plot_circAD

        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')
        t.show_text('text' => '\textbf{Z = 0.02}', 'color' => Black, 'x' => 0.85, 'y' => 0.1)
	background

        xdata1 = galr_6.safe_log10
	#xdata2 = (rhowimp_6/v_star_6).safe_log10
	ydata = (lwimp0_6/(lpp0_6+lCNO0_6)).safe_log10
	xdata3 = galr_8.safe_log10
	ydata3 = (lwimp0_8/(lpp0_8+lCNO0_8)).safe_log10
	aryx = [xdata1, xdata3]
	aryy = [ydata, ydata3]

	#    t.subplot do
	#	t.xaxis_loc = t.xlabel_side = TOP;
	#	t.left_edge_type = AXIS_HIDDEN        	
	#	t.bottom_edge_type = AXIS_HIDDEN
	#	t.right_edge_type = AXIS_HIDDEN		
	#	xmax = xdata2.max; xmin = xdata2.min; xmargin = @margin*(xmax-xmin)
        #	right = xmax + xmargin; left = xmin - xmargin 
        #	ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        #	top = ymax + ymargin; bottom = ymin - ymargin
	#	t.show_xlabel('$\frac{\rho_\chi}{v_\star}$ ($\log_{10}$ GeV s km$^{-1}$)')
        #	t.show_plot('left_boundary' => left, 'right_boundary' => right, 'top_boundary' => top, 'bottom_boundary' => bottom) do end
	#    end

        xmax = max_of_many(aryx); xmin = min_of_many(aryx); xmargin = @margin*(xmax-xmin)
        right = xmax + xmargin; left = xmin - xmargin 
        ymax = max_of_many(aryy); ymin = min_of_many(aryy); ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
        t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.05,'plot_right_margin' => 0,
         'legend_top_margin' => 0.02) do
	    t.subplot do
	#	t.top_edge_type = AXIS_HIDDEN	
		t.show_plot('left_boundary' => left, 'right_boundary' => right, 'top_boundary' => top, 'bottom_boundary' => bottom) do
			t.show_xlabel('Orbital radius, $\log_{10}\big(\frac{r}{\mathrm{pc}}\big)$')	
			colours.sort.reverse.each do |m, col| 
                		sorted_data = [Array.new, Array.new, Array.new]
                		sorted_data[0] = (0..mass_6.length-1).find_all{|i| (mass_6[i]/m-1).abs < Epsilon}
                		sorted_data[0].each_index do |i|
                			sorted_data[1][i] = xdata1[sorted_data[0][i]]
                			sorted_data[2][i] = ydata[sorted_data[0][i]]
                			#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
				end
				if sorted_data[1][0] != nil then
					x = sorted_data[1].sort
					y = sorted_data[2].sort {|a,b| sorted_data[1][sorted_data[2].index(a)] <=> sorted_data[1][sorted_data[2].index(b)]}
					t.show_polyline(x, y, col, m.to_s+'\,M$_\odot$')
				end
			end
			colours.sort.reverse.each do |m, col| 
        	        	sorted_data = [Array.new, Array.new, Array.new]
        	        	sorted_data[0] = (0..mass_8.length-1).find_all{|i| (mass_8[i]/m-1).abs < Epsilon}
        	        	sorted_data[0].each_index do |i|
        	        		sorted_data[1][i] = xdata3[sorted_data[0][i]]
        	        		sorted_data[2][i] = ydata3[sorted_data[0][i]]
        	        	end
				if sorted_data[1][0] != nil then
					x = sorted_data[1].sort
					y = sorted_data[2].sort {|a,b| sorted_data[1][sorted_data[2].index(a)] <=> sorted_data[1][sorted_data[2].index(b)]}
					t.show_polyline(x, y, col, m.to_s+'\,M$_\odot$, $f(u)$ truncated', Line_Type_Dash) 
				end
	       	        	t.line_type = Line_Type_Solid           	
			end
		end
	    end
        end       
    end

    def plot_circNFW

        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')
        t.show_text('text' => '\textbf{Z = 0.02}', 'color' => Black, 'x' => 0.85, 'y' => 0.1)
	background	        
	
        xdata1 = galr_7.safe_log10
	xdata2 = (rhowimp_7/v_star_7).safe_log10
	ydata = (lwimp0_7/(lpp0_7+lCNO0_7)).safe_log10
	#    t.subplot do
	#	t.xaxis_loc = t.xlabel_side = TOP;
	#	t.left_edge_type = AXIS_HIDDEN        	
	#	t.bottom_edge_type = AXIS_HIDDEN
	#	t.right_edge_type = AXIS_HIDDEN		
	#	xmax = xdata2.max; xmin = xdata2.min; xmargin = @margin*(xmax-xmin)
        #	right = xmax + xmargin; left = xmin - xmargin 
        #	ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        #	top = ymax + ymargin; bottom = ymin - ymargin
	#	t.show_xlabel('$\frac{\rho_\chi}{v_\star}$ ($\log_{10}$ GeV s km$^{-1}$)')
        #	t.show_plot('left_boundary' => left, 'right_boundary' => right, 'top_boundary' => top, 'bottom_boundary' => bottom) do end
	#    end

        xmax = xdata1.max; xmin = xdata1.min; xmargin = @margin*(xmax-xmin)
        right = xmax + xmargin; left = xmin - xmargin 
        ymax = ydata.max; ymin = ydata.min; ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
        t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.05,'plot_right_margin' => 0,
            'legend_top_margin' => 0.02) do
	    t.subplot do
	#	t.top_edge_type = AXIS_HIDDEN	
		t.show_plot('left_boundary' => left, 'right_boundary' => right, 'top_boundary' => top, 'bottom_boundary' => bottom) do
			t.show_xlabel('Orbital radius, $\log_{10}\big(\frac{r}{\mathrm{pc}}\big)$')	
			colours.sort.reverse.each do |m, col| 
                		sorted_data = [Array.new, Array.new, Array.new]
                		sorted_data[0] = (0..mass_7.length-1).find_all{|i| (mass_7[i]/m-1).abs < Epsilon}
                		sorted_data[0].each_index do |i|
                			sorted_data[1][i] = xdata1[sorted_data[0][i]]
                			sorted_data[2][i] = ydata[sorted_data[0][i]]
                			#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
				end
				if sorted_data[1][0] != nil then
					x = sorted_data[1].sort
					y = sorted_data[2].sort {|a,b| sorted_data[1][sorted_data[2].index(a)] <=> sorted_data[1][sorted_data[2].index(b)]}
					t.show_polyline(x, y, col, m.to_s+'\,M$_\odot$')
				end
                	end
		end
	    end
        end       

    end


    def plot_elliptical_10yr
	plot_elliptical_general(10)
    end

    def plot_elliptical_50yr
	plot_elliptical_general(50)
    end

    def plot_elliptical_cpc
	plot_elliptical_general(0.01)
    end


    def plot_elliptical_general(orb)

        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('Mean initial capture rate, $\log_{10}\big(\frac{\langle C\rangle_{5\,\mathrm{orbits}}}{n\,\mathrm{yr}^{-1}}\big)$')
	if (orb == 0.01) then
		msg = '$\mathbf{r_{max}=0.01}$\,\textbf{pc}'
	else
		msg = '$\mathbf{P='+orb.to_s+'}$\,\textbf{yr}'
	end
        t.show_text('text' => msg, 'color' => Black, 'x' => 0.75, 'y' => 0.2)
	t.show_text('text' => '\textbf{Z = 0.02}', 'color' => Black, 'x' => 0.75, 'y' => 0.1)
	background	        
	
        xdata1 = (1-e_9).safe_log10
	xdata2 = (1-e_10).safe_log10
	ydata1 = cap_9.safe_log10
	ydata2 = cap_10.safe_log10
	aryx = [xdata1, xdata2]
	aryy = [ydata1, ydata2]

        xmax = max_of_many(aryx); xmin = min_of_many(aryx); xmargin = @margin*(xmax-xmin)
        #right = xmax + xmargin; left = xmin - xmargin 
	ymax = max_of_many(aryy); ymin = min_of_many(aryy); ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
        right = 0.25; left = -3.5; top = 44;       
	t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.05,'plot_right_margin' => 0,
            'legend_top_margin' => 0.42) do
		t.show_plot('left_boundary' => left, 'right_boundary' => right, 'top_boundary' => top, 'bottom_boundary' => bottom) do
			t.show_xlabel('Modified orbital ellipticity, $\log_{10}(1-e)$')
			colours.sort.reverse.each do |m, col| 
                		sorted_data = [Array.new, Array.new, Array.new]
                		sorted_data[0] = (0..mass_9.length-1).find_all{|i| (mass_9[i]/m-1).abs < Epsilon && orbit_9[i] == orb}
                		sorted_data[0].each_index do |i|
                			sorted_data[1][i] = xdata1[sorted_data[0][i]]
                			sorted_data[2][i] = ydata1[sorted_data[0][i]]
                			#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
				end
				if sorted_data[1][0] != nil then
					x = sorted_data[1].sort
					y = sorted_data[2].sort {|a,b| sorted_data[1][sorted_data[2].index(a)] <=> sorted_data[1][sorted_data[2].index(b)]}
					t.show_polyline(x, y, col, m.to_s+'\,M$_\odot$, AC+spike', Line_Type_Solid)
				end
			end
			colours.sort.reverse.each do |m, col| 
				sorted_data = [Array.new, Array.new, Array.new]
                		sorted_data[0] = (0..mass_10.length-1).find_all{|i| (mass_10[i]/m-1).abs < Epsilon && orbit_10[i] == orb}
                		sorted_data[0].each_index do |i|
                			sorted_data[1][i] = xdata2[sorted_data[0][i]]
                			sorted_data[2][i] = ydata2[sorted_data[0][i]]
                			#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
				end
				if sorted_data[1][0] != nil then
					x = sorted_data[1].sort
					y = sorted_data[2].sort {|a,b| sorted_data[1][sorted_data[2].index(a)] <=> sorted_data[1][sorted_data[2].index(b)]}
					t.show_polyline(x, y, col, m.to_s+'\,M$_\odot$, NFW+spike', Line_Type_Dash)
				end				
			end
		end
	end       

    end


 def plot_elliptical_halo

	m = 1.0
	orb = 10
	col1 = Burgundy
        col2 = Goldenrod
        col3 = SeaGreen
        col4 = BrightBlue
        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
        t.show_ylabel('Mean initial capture rate, $\log_{10}\big(\frac{\langle C\rangle_{5\,\mathrm{orbits}}}{n\,\mathrm{yr}^{-1}}\big)$')
	msg = '$\mathbf{P=10}$\,\textbf{yr}'
	t.show_text('text' => msg, 'color' => Black, 'x' => 0.25, 'y' => 0.2)
	t.show_text('text' => '\textbf{Z = 0.02}', 'color' => Black, 'x' => 0.25, 'y' => 0.1)
	background	        
	
        xdata1 = (1-e_9).safe_log10
	xdata2 = (1-e_13).safe_log10
	xdata3 = (1-e_11).safe_log10
	xdata4 = (1-e_12).safe_log10
	ydata1 = cap_9.safe_log10
	ydata2 = cap_13.safe_log10
	ydata3 = cap_11.safe_log10
	ydata4 = cap_12.safe_log10
	aryx = [xdata1, xdata2, xdata3, xdata4]
	aryy = [ydata1, ydata2, ydata3, ydata4]

        xmax = max_of_many(aryx); xmin = min_of_many(aryx); xmargin = @margin*(xmax-xmin)
        #right = xmax + xmargin; left = xmin - xmargin 
	ymax = max_of_many(aryy); ymin = min_of_many(aryy); ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
        right = 0.25; left = -3.5; top = 44; bottom = 35       
	t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.05,'plot_right_margin' => 0,
            'legend_top_margin' => 0.38) do
		t.show_plot('left_boundary' => left, 'right_boundary' => right, 'top_boundary' => top, 'bottom_boundary' => bottom) do
			t.show_xlabel('Modified orbital ellipticity, $\log_{10}(1-e)$')

			sorted_data = [Array.new, Array.new, Array.new]
                	sorted_data[0] = (0..mass_12.length-1).find_all{|i| (mass_12[i]/m-1).abs < Epsilon && orbit_12[i] == orb}
                	sorted_data[0].each_index do |i|
                		sorted_data[1][i] = xdata4[sorted_data[0][i]]
                		sorted_data[2][i] = ydata4[sorted_data[0][i]]
                		#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
			end
			if sorted_data[1][0] != nil then
				x = sorted_data[1].sort
				y = sorted_data[2].sort {|a,b| sorted_data[1][sorted_data[2].index(a)] <=> sorted_data[1][sorted_data[2].index(b)]}
				t.show_polyline(x, y, col1, m.to_s+'\,M$_\odot$, AC+spike, full $N$-body', Line_Type_Solid)
			end				

			sorted_data = [Array.new, Array.new, Array.new]
                	sorted_data[0] = (0..mass_11.length-1).find_all{|i| (mass_11[i]/m-1).abs < Epsilon && orbit_11[i] == orb}
                	sorted_data[0].each_index do |i|
                		sorted_data[1][i] = xdata3[sorted_data[0][i]]
                		sorted_data[2][i] = ydata3[sorted_data[0][i]]
                		#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
			end
			if sorted_data[1][0] != nil then
				x = sorted_data[1].sort
				y = sorted_data[2].sort {|a,b| sorted_data[1][sorted_data[2].index(a)] <=> sorted_data[1][sorted_data[2].index(b)]}
				t.show_polyline(x, y, col2, m.to_s+'\,M$_\odot$, AC+spike, truncated $N$-body', Line_Type_Dot)
			end

			sorted_data = [Array.new, Array.new, Array.new]
                	sorted_data[0] = (0..mass_13.length-1).find_all{|i| (mass_13[i]/m-1).abs < Epsilon && orbit_13[i] == orb}
                	sorted_data[0].each_index do |i|
                		sorted_data[1][i] = xdata2[sorted_data[0][i]]
                		sorted_data[2][i] = ydata2[sorted_data[0][i]]
                		#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
			end
			if sorted_data[1][0] != nil then
				x = sorted_data[1].sort
				y = sorted_data[2].sort {|a,b| sorted_data[1][sorted_data[2].index(a)] <=> sorted_data[1][sorted_data[2].index(b)]}
				t.show_polyline(x, y, col3, m.to_s+'\,M$_\odot$, AC+spike, full isothermal', Line_Type_Solid)
			end

			sorted_data = [Array.new, Array.new, Array.new]
                	sorted_data[0] = (0..mass_9.length-1).find_all{|i| (mass_9[i]/m-1).abs < Epsilon && orbit_9[i] == orb}
                	sorted_data[0].each_index do |i|
                		sorted_data[1][i] = xdata1[sorted_data[0][i]]
                		sorted_data[2][i] = ydata1[sorted_data[0][i]]
                		#mark_spot(sorted_data[1][i], sorted_data[2][i], colours[sorted_data[0][i]])
			end
			if sorted_data[1][0] != nil then
				x = sorted_data[1].sort
				y = sorted_data[2].sort {|a,b| sorted_data[1][sorted_data[2].index(a)] <=> sorted_data[1][sorted_data[2].index(b)]}
				t.show_polyline(x, y, col4, m.to_s+'\,M$_\odot$, AC+spike, truncated isothermal', Line_Type_Dot)
			end
		end
	end       

    end


def plot_elliptical_lwimp

	ellipticities = [0.0, -0.301030, -1.00000, -2.00000, -3.00000, -3.30103]
	lwAC06 = [-3.3089310, -1.1148967, 0.50869499, 0.72394085, 1.0270967, 1.0373444]
	lwAC10 = [-3.5332096, -1.8308621, -0.15877802, 1.5181346]
	ellipticities_lwAC10 = [0.0, -0.301030, -1.00000, -2.00000]
	lwAC15 = [-3.7127568, -2.2888648, -0.43732893, 1.0450190, 1.7263796, 1.8612660]
	lwNFW06 = [-5.1668241,-2.9727899, -1.3852531, -0.79437481, -0.40280933, -0.42072470]
	lwNFW10 = [-5.3433968, -3.6410492, -1.9593876, -0.57027710, -0.21158129, -0.18799070]
	lwNFW15 = [-5.5243786, -4.1004866, -2.2447689, -0.99597213, -0.64357206, -0.62828277]

        t.legend_text_dy = 1.3
        t.ylabel_scale = 1.2
        t.xlabel_scale = 1.2
        t.yaxis_numeric_label_scale = 1
        t.xaxis_numeric_label_scale = 1
	t.show_ylabel('WIMP luminosity, $\log_{10}\Big[\frac{L_\mathrm{W,max}}{L_\mathrm{nuc}(0)}\Big]$')
        msg = '$\mathbf{P=10}$\,\textbf{yr}'
	t.show_text('text' => msg, 'color' => Black, 'x' => 0.75, 'y' => 0.2)
	t.show_text('text' => '\textbf{Z = 0.02}', 'color' => Black, 'x' => 0.75, 'y' => 0.1)
	background	        
	
	aryy = [lwAC06, lwAC10, lwAC15, lwNFW06, lwNFW06, lwNFW15]

        xmax = max_of_many(ellipticities); xmin = min_of_many(ellipticities); xmargin = @margin*(xmax-xmin)
        ymax = max_of_many(aryy); ymin = min_of_many(aryy); ymargin = @margin*(ymax-ymin)
        top = ymax + ymargin; bottom = ymin - ymargin
        right = 0.25; left = -3.5; top = 3; bottom = -7      
	t.show_plot_with_legend('legend_scale' => 1.4,'legend_left_margin' => 0.05,'plot_right_margin' => 0,
            'legend_top_margin' => 0.5) do
		t.show_plot('left_boundary' => left, 'right_boundary' => right, 'top_boundary' => top, 'bottom_boundary' => bottom) do
			t.show_xlabel('Modified orbital ellipticity, $\log_{10}(1-e)$')
			t.show_polyline(ellipticities, lwAC06, colours[0.6], '0.6\,M$_\odot$, AC+spike', Line_Type_Solid)
			t.show_polyline(ellipticities_lwAC10[0,3], lwAC10[0,3], colours[1.0], '1.0\,M$_\odot$, AC+spike', Line_Type_Solid)
			t.show_polyline(ellipticities, lwAC15, colours[1.5], '1.5\,M$_\odot$, AC+spike', Line_Type_Solid)
			t.show_polyline(ellipticities, lwNFW06, colours[0.6], '0.6\,M$_\odot$, NFW+spike', Line_Type_Dash)
			t.show_polyline(ellipticities, lwNFW10, colours[1.0], '1.0\,M$_\odot$, NFW+spike', Line_Type_Dash)
			t.show_polyline(ellipticities, lwNFW15, colours[1.5], '1.5\,M$_\odot$, NFW+spike', Line_Type_Dash)
		t.show_arrow('head' => [ellipticities[3],lwAC10[3]], 'tail'=> [ellipticities[2],lwAC10[2]], 'head_scale' => 1.0, 'head_color' => colours[1.0], 'tail_marker' => 'None', 'line_color' => colours[1.0])
		end
	end
end


    def add_info_line(xs, ys)
        t.line_type = Line_Type_Dash
        t.line_color = SlateGray
        t.line_width = 1
        t.append_points_to_path(xs, ys)
        t.stroke
    end

    def min_of_many(ary, y_limit = nil)
        return nil if ary == nil || ary.size == 0
        ymin = Dvector.min_of_many(ary)
        ymin = y_limit if y_limit != nil && ymin < y_limit
        return ymin
    end

    def max_of_many(ary, y_limit = nil)
        return nil if ary == nil || ary.size == 0
        ymax = Dvector.max_of_many(ary)
        ymax = y_limit if y_limit != nil && ymax > y_limit
        return ymax
    end

    def mark_spot(x, y, color = Red)
        t.show_marker('x' => x, 'y' => y, 'marker' => Bullet, 'scale' => 0.5, 'color' => color);
    end

    def read_ZAMS(path = nil)
        path = '' if path == nil
        path = path + '/' if path.length > 0 && path[-1..-1] != '/'
        @zams = [
            @zams_log_Center_Density = Dvector.new,
            @zams_log_Center_Temp = Dvector.new,
            @zams_log_Luminosity = Dvector.new,
            @zams_log_Surface_Temp = Dvector.new,
            @zams_Ms = Dvector.new ]
        Dvector.read(path + 'ZAMS.data', @zams, 2)
    end

    def stroke_track(xs, ys, color=Black, legend=nil, type=Line_Type_Solid, track_first=@track_first)
        t.stroke_color = color
        t.line_type = type
        if @gap_in_tracks
            t.append_points_with_gaps_to_path(xs,ys,[@helium_gap_start-track_first],false)
        else
            t.append_points_to_path(xs,ys)
        end
        t.save_legend_info(legend) if legend != nil
        t.stroke
    end

end

