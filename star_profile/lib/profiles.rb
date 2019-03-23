# profiles.rb

class ProfilePlots

    include Math
    include Tioga
    include FigureConstants
    
    def t
        @figure_maker
    end

    def d
        @profile_data
    end
    
    def set_profile_data(data)
        @profile_data = data
        @have_data = true
        t.need_to_reload_data = false
    end

    def hr
        @hr_trho_data
    end

    def eos
        @eos_data
    end
    
    def speeds
        read_data
        show_model_number
        xs = d.sx_logP
        t.rescale(0.8)
        t.do_box_labels('Ratio of Convection Velocity to Sound Speed', 'logP', 'convection/sound')
        ys = d.sx_CV  / d.sx_CS
        ys.safe_log10!
        t.show_plot('left_boundary' => xs.max, 'right_boundary' => xs.min-1,
            'top_boundary' => 0.2, 'bottom_boundary' => -5) { t.show_polyline(xs, ys, Blue) }
    end
    
    def initialize(logfile, savedir = 'profiles_out', slave_mode = false)
    
        @figure_maker = FigureMaker.default
        @rescale_value = 1
        t.save_dir = savedir unless savedir == nil
        t.def_eval_function { |str| eval(str) }
        t.tex_preview_preamble = t.tex_preview_preamble + "\n\t\\include{color_names}\n"
        
        if not slave_mode
            @profile_name = logfile
            #t.def_figure("speeds") { speeds }
            t.def_figure("Abundances_by_Mass") { abundances_by_mass }
            t.def_figure("Power_by_Mass") { power_by_mass }
            t.def_figure("Ratios_by_Mass") { ratios_by_mass }
            t.def_figure("Trio_by_Mass") { trio_by_mass }
            t.def_figure("Convection_by_logP") { convection_by_logP }
            t.def_figure("Ionization_by_logP") { ionization_by_logP }
            t.def_figure("Ratios_by_logP") { ratios_by_logP }
            t.def_figure("Trio_by_logP") { trio_by_logP }

            #t.def_figure("Duo_by_Mass") { duo_by_mass }
            t.def_figure("Convection_by_Mass") { convection_by_mass }
            #t.def_figure("Superadiabaticity_by_Mass") { grad_star_by_mass } # show where super-adiabatic
            t.def_figure("Grad_star_by_logP") { grad_star_by_logP }
            t.def_figure("Abundances_by_logP") { abundances_by_logP }
            #t.def_figure("Power_by_logP") { power_by_logP }
            #t.def_figure("Trio2_by_logP") { trio2_by_logP }
            #t.def_figure("Abundances_Both_Ways") { abundances_by_both }
            #t.def_figure("Convection_Both_Ways") { convection_by_both }
            #t.def_figure("Trio_Both_Ways") { trio_by_both }
            
            t.def_figure("Full_Profile") { full_profile }

            @energy_eos_plot = t.def_figure("Energy_EoS") { energy_eos }
            @pressure_eos_plot = t.def_figure("Pressure_EoS") { pressure_eos }
            @opacity_eos_plot = t.def_figure("Opacity_EoS") { opacity_eos }
            @ionizaton_eos_plot = t.def_figure("Ionization_EoS") { ionization_eos }
            @entropy_eos_plot = t.def_figure("Entropy_EoS") { entropy_eos }
            #@psi_eos_plot = t.def_figure("PSI_EoS") { psi_eos }
            
            t.def_figure("Entropy_and_History") { entropy_eos_and_history }
            t.def_figure("Trio2_and_History") { trio2_and_history }
            t.def_figure("Trio2_by_Mass") { trio2_by_mass }
            #t.def_figure("HR_history") { plot_H_R_T_RHO }
        
            @have_hr_data = false
    
            read_PSIs('star_data')
            read_ZAMS('star_data')
            @hr_name = "../run/HR_history.log"
            @hr_trho_data = HR_History_Data.new
        
            @profile_data = ProfileData.new
            t.auto_refresh_filename = @profile_name
        
        end

        @eos_data = EosData.new

        @mass_50_color = LightSteelBlue
        @mass_95_color = MediumSlateBlue
        @mass_999_color = MediumBlue
        
        @show_mesh_flag = true

        @have_data = false
        @mass_xlabel = 'enclosed mass \small{($\mathrm{M_{\odot}} $)}'
        @log_mass_xlabel = 'Enclosed Mass \small{($log$ $\mathrm{M_{\odot}} $)}'
        @logP_xlabel = '$ \log P $'
        
        @subplot_scale = 0.8
        @subplot_margin = 0.01
        @half_plot = 0.38  
        @legend_scale = 0.95
        @legend_on_right = true
        @plot_right = 0.15
        @legend_left = 0.7
        @legend_top = 0.2
        @trio_margin = 0.02
        @side_by_side_margin = 0.51
        @plot_with_legend_right_margin = 0.10
        @legend_left_margin = 0.92
        @legend_scale = 0.9
    
        @plot_with_legend_dict = { 'plot_right' => @plot_right, 'legend_left' => @legend_left, 'legend_top' => @legend_top }

        @rho_T_plot = lambda {|xs, title, xlabel, xleft, xright| rho_T_2Ys(xs, title, xlabel, xleft, xright) }
        @density_plot = lambda {|xs, title, xlabel, xleft, xright| density(xs, title, xlabel, xleft, xright) }
        @temperature_plot = lambda {|xs, title, xlabel, xleft, xright| temperature(xs, title, xlabel, xleft, xright) }
        @abundances_plot = lambda {|xs, title, xlabel, xleft, xright| abundances(xs, title, xlabel, xleft, xright) }
        @convection_plot = lambda {|xs, title, xlabel, xleft, xright| convection(xs, title, xlabel, xleft, xright) }
        @ionization_plot = lambda {|xs, title, xlabel, xleft, xright| ionization(xs, title, xlabel, xleft, xright) }
        @ratios_plot = lambda {|xs, title, xlabel, xleft, xright| ratios(xs, title, xlabel, xleft, xright) }
        @power_plot = lambda {|xs, title, xlabel, xleft, xright| power(xs, title, xlabel, xleft, xright) }
        @grad_star_plot = lambda {|xs, title, xlabel, xleft, xright| grad_star(xs, title, xlabel, xleft, xright) }

        @temp_vec1 = Dvector.new
        @temp_vec2 = Dvector.new
        
        t.rescale(@rescale_value)
        t.def_enter_page_function { enter_page }
            
    end
    
    def enter_page
        t.page_setup(@rescale_value*11*72/2,@rescale_value*8.5*72/2)
        t.set_frame_sides(0.15,0.85,0.85,0.15) # left, right, top, bottom in page coords        
    end
    
    def read_data
        if !t.root_figure || (@have_data && !t.need_to_reload_data)
            return
        end
        d.read_profile(@profile_name)
        @have_data = true
        t.model_number = d.model_Number
        t.need_to_reload_data = false
    end
    
end
