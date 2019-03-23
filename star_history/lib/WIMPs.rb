# WIMPs.rb

class StarHistory

    include Math
    include FigureConstants
    
    def get_track(xs)
        return xs[@track_first .. @track_last]
    end

    def special_boundaries(xs,ys,margin,ymin=nil,ymax=nil,ymargin=nil)
        extra_on_top = 2.5
        extra_on_bottom = 2.5
        ymin = ys.min if ymin == nil
        ymax = ys.max if ymax == nil
        height = (ymax == ymin)? 1 : ymax - ymin
        ymargin = margin if ymargin == nil
        top_boundary = ymax + extra_on_top * ymargin * height
        bottom_boundary = ymin - extra_on_bottom * ymargin * height
        return [ @start_age, @end_age, top_boundary, bottom_boundary ]
    end

    def plot_one_value_WIMPs(ys, ylabel1, ylabel2 = '', ymin_limit = nil, legend_size = false)
        xs = d.star_Age[@track_first .. @track_last]
        ys = ys[@track_first .. @track_last]
        t.yaxis_numeric_label_scale = 0.5
	t.yaxis_major_tick_length = 0.4
        t.yaxis_minor_tick_length = 0.2
	t.yaxis_major_tick_width = 0.5
	t.yaxis_minor_tick_width = 0.5
	t.xaxis_major_tick_length = 0.4
        t.xaxis_minor_tick_length = 0.2
	t.xaxis_major_tick_width = 0.5
	t.xaxis_minor_tick_width = 0.5
	t.yaxis_numeric_label_scale = 0.75
        t.xaxis_numeric_label_scale = 0.8
        t.xlabel_scale = 1.2
        ymin = ys.min
        ymin = ymin_limit if ymin_limit != nil && ymin < ymin_limit
        if legend_size then
            t.show_text('text' => ylabel1, 'color' => Black, 'x' => 0.28, 'y' => 0.6, 'scale' => 1.0)	
	    t.show_text('text' => ylabel2, 'color' => Black, 'x' => 0.28, 'y' => 0.35, 'scale' => 1.0)	
            t.show_plot_with_legend('plot_left_margin' => 0.4) do
                t.show_plot('boundaries' => special_boundaries(xs,ys,@margin,ymin)) do
                    #background
                    #mark_profiles_on_x
                    stroke_track(xs,ys,Black)#DarkOrange
                    t.show_xlabel(@age_xlabel)
                end
            end
        else
            t.show_ylabel(ylabel); t.no_ylabel
            t.show_plot('boundaries' => plot_boundaries(xs,ys,@margin,ymin)) do
                background
                mark_profiles_on_x
                stroke_track(xs,ys,DarkOrange)
                mark_line_at_loc(xs, ys, -1, color = BrightBlue)
            end
        end
    end

    def plot_luminosities(xs)
	t.yaxis_major_tick_length = 0.4
        t.yaxis_minor_tick_length = 0.2
	t.yaxis_major_tick_width = 0.5
	t.yaxis_minor_tick_width = 0.5
	t.xaxis_major_tick_length = 0.4
        t.xaxis_minor_tick_length = 0.2
	t.xaxis_major_tick_width = 0.5
	t.xaxis_minor_tick_width = 0.5       
        t.legend_text_dy = 1.3
	#t.show_ylabel('Power ($\log_{10}$ $L/\mathrm{L}_\odot$)'); t.no_ylabel
        t.show_text('text' => 'Power', 'color' => Black, 'x' => 0.28, 'y' => 0.6, 'scale' => 1.0)	
	t.show_text('text' => '$\log_{10}(L/\mathrm{L}_\odot)$', 'color' => Black, 'x' => 0.28, 'y' => 0.35, 'scale' => 1.0)	
	sun_L = 3.844e33      
        h = get_track_log(d.power_H_burn)
        other = get_track_log(d.power_He_burn + d.power_Metal_burn)
        wimps = get_track_log(d.wimp_Energy_Inject/sun_L)
        total = get_track_log(d.power_H_burn + d.power_He_burn + d.power_Metal_burn + d.wimp_Energy_Inject/sun_L)
        t.yaxis_numeric_label_scale = 0.75
        t.xaxis_numeric_label_scale = 0.8
        min_val = min_of_many([(d.wimp_Energy_Inject[2..-1]/sun_L).safe_log10.min, d.power_H_burn.safe_log10.min])
        #t.show_plot_with_legend('legend_left_margin' => 0.62,'legend_top_margin' => 0.23, 'plot_left_margin' => 0.4) do
        t.subplot('left_margin' => 0.4, 'right_margin' => 0.18) do
           t.show_plot('boundaries' => get_boundaries([ h, other, total, wimps ], min_val)) do
                #background
                #mark_profiles_on_x
                stroke_track(xs, wimps, Black, 'WIMPs')
                stroke_track(xs, h, Goldenrod, 'H (pp $+$ CNO)')
                #stroke_track(xs, other, SkyBlue, 'He $+$ metals')
                stroke_track(xs, total, FireBrick, 'Total burning', Line_Type_Dot)
                t.show_xlabel(@age_xlabel)
            end
        end
    end

    def plot_Tc2
        plot_one_value_WIMPs(d.log_center_Temp, '$\log_{10}$ $T_\mathrm{c}$', '', nil, true)
    end
    
    def plot_Rho2
        plot_one_value_WIMPs(d.log_center_Density, '$\log_{10}$ $\rho_\mathrm{c}$', '', nil, true)
    end

    def plot_Pc        
        plot_one_value_WIMPs(d.log_center_Pressure, '$\log_{10}$ $P_\mathrm{c}$', '', nil, true)
    end

    def plot_n_WIMPs        
        plot_one_value_WIMPs(d.n_WIMPs.safe_log10, 'Total WIMPs', '$\log_{10}(n)$', nil, true)#, d.n_WIMPs[1].safe_log10)
    end

    def plot_cap
        plot_one_value_WIMPs(d.cap.safe_log10, 'Capture rate', '$\log_{10}\big(\frac{C}{n\,\mathrm{yr}^{-1}}\big)$', nil, true)
    end

    def plot_ann
        plot_one_value_WIMPs(d.ann_Times2.safe_log10, 'Annihilation rate', '$\log_{10}\big(\frac{A}{n\,\mathrm{yr}^{-1}}\big)$', d.cap.safe_log10.min, true)
    end

    def plot_dNdt
        plot_one_value_WIMPs(d.dNdt.abs.safe_log10, 'Net change', '$\log_{10}\big(\frac{|\Delta n|}{n\,\mathrm{yr}^{-1}}\big)$', nil, true)
    end

    def plot_r_Gal
        plot_one_value_WIMPs(d.r_Gal.safe_log10, 'Orbital distance', '$\log_{10}\big(\frac{r}{\mathrm{pc}}\big)$', nil, true)
    end

    def plot_rho_WIMP
        plot_one_value_WIMPs(d.rho_WIMP.safe_log10, 'WIMP density', '$\log_{10}\big(\frac{\rho_\chi}{\mathrm{GeV\,cm}^{-3}}\big)$', nil, true)
    end

    def plot_v_Star
        plot_one_value_WIMPs(d.v_Star.safe_log10, 'Stellar velocity', '$\log_{10}\big(\frac{v_\star}{\mathrm{km\,s}^{-1}}\big)$', nil, true)
    end

    def plot_knudsen(xs)
        t.show_ylabel('Knudsen Parameter'); t.no_ylabel
        k = get_track(d.knudsen)
        t.show_plot_with_legend do
            t.show_plot('boundaries' => get_boundaries([ k ], 0)) do
                background
                mark_profiles_on_x
                stroke_track(xs, k, Black)
                t.show_xlabel(@age_xlabel)
            end
        end
    end

    def plot_r_chi(xs)
        t.show_ylabel('$\log_{10}$ normalised $r_\mathrm{\chi}$'); t.no_ylabel
        r_Sun = 0.69598e11
        r1 = get_track(d.r_chi.safe_log10 - d.log_Radius - r_Sun.safe_log10)
        r2 = get_track_log(d.r_chi/r_Sun)
        t.show_plot_with_legend do
            t.show_plot('boundaries' => get_boundaries([ r1, r2 ], -20)) do
                background
                mark_profiles_on_x
                stroke_track(xs, r1, BrightBlue, '$r_\chi$ / $\mathrm{R}_\star$')
                stroke_track(xs, r2, SkyBlue, '$r_\chi$ / $\mathrm{R}_\odot$')
                t.show_xlabel(@age_xlabel)
            end
        end
    end

    def plot_temps(xs)
        t.show_ylabel('$\log_{10}$ $T$ (K)'); t.no_ylabel
        tc = get_track(d.log_center_Temp)
        tw = get_track_log(d.Tw)
        t.show_plot_with_legend do
            t.show_plot('boundaries' => get_boundaries([ tc, tw ], 0)) do
                background
                mark_profiles_on_x
                stroke_track(xs, tc, Goldenrod, '$T_\mathrm{c}$')
                stroke_track(xs, tw, FireBrick, '$T_\mathrm{W}$')
                t.show_xlabel(@age_xlabel)
            end
        end
    end

    def plot_densities(xs)
        t.show_ylabel('WIMP density ($\log_{10}$ GeV/cm$^{3}$)'); t.no_ylabel
        cen = get_track_log(d.centre_WIMP_Dens)
        rchi = get_track_log(d.r_chi_WIMP_Dens)
        t.show_plot_with_legend do
            t.show_plot('boundaries' => get_boundaries([ cen, rchi ], -5)) do
                background
                mark_profiles_on_x
                stroke_track(xs, cen, BrightBlue, '$r=0$')
                stroke_track(xs, rchi, Goldenrod, '$r=r_\chi$')
                t.show_xlabel(@age_xlabel)
            end
        end
    end
    
    
    def plot_WIMP_luminosity

        setup_data
        set_track_for_plots
        t.legend_text_dy = 1.3
        t.rescale(0.8)
        t.ylabel_scale = 0.8
        t.yaxis_numeric_label_scale = 0.6
        num_plots = 4; row = 1
        xs = d.star_Age[@track_first .. @track_last]
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 1)) do # first
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.no_xlabel
            plot_luminosities(xs)
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            t.no_xlabel
            plot_Tc2
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            t.no_xlabel
            plot_Rho2
        row += 1
        end
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.top_edge_type = AXIS_HIDDEN
            plot_Pc
        end
    end


    def plot_WIMP_combo

        setup_data
        set_track_for_plots
        t.legend_text_dy = 1
        t.legend_scale = 0.5
	t.rescale(0.8)
        t.ylabel_scale = 0.7
        t.yaxis_numeric_label_scale = 0.5
        t.xlabel_scale = 0.7
        t.xaxis_numeric_label_scale = 0.5
        num_plots = 7; row = 1
        xs = d.star_Age[@track_first .. @track_last]
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 1)) do # first
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.yaxis_type = AXIS_WITH_MAJOR_TICKS_AND_NUMERIC_LABELS
            t.right_edge_type = AXIS_WITH_MAJOR_TICKS_ONLY
            t.top_edge_type = AXIS_LINE_ONLY
            t.no_xlabel
            plot_luminosities(xs)
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do # first
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.yaxis_type = AXIS_WITH_MAJOR_TICKS_ONLY
            t.right_edge_type = AXIS_WITH_MAJOR_TICKS_AND_NUMERIC_LABELS
            t.top_edge_type = AXIS_HIDDEN
            t.no_xlabel
            plot_n_WIMPs
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.yaxis_type = AXIS_WITH_MAJOR_TICKS_AND_NUMERIC_LABELS
            t.right_edge_type = AXIS_WITH_MAJOR_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            t.no_xlabel
            plot_cap
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.yaxis_type = AXIS_WITH_MAJOR_TICKS_ONLY
            t.right_edge_type = AXIS_WITH_MAJOR_TICKS_AND_NUMERIC_LABELS
            t.top_edge_type = AXIS_HIDDEN
            t.no_xlabel
            plot_ann
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.yaxis_type = AXIS_WITH_MAJOR_TICKS_AND_NUMERIC_LABELS
            t.right_edge_type = AXIS_WITH_MAJOR_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            t.no_xlabel
            plot_r_Gal
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.yaxis_type = AXIS_WITH_MAJOR_TICKS_ONLY
            t.right_edge_type = AXIS_WITH_MAJOR_TICKS_AND_NUMERIC_LABELS
            t.top_edge_type = AXIS_HIDDEN
            t.no_xlabel
            plot_rho_WIMP
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.top_edge_type = AXIS_HIDDEN
            t.yaxis_type = AXIS_WITH_MAJOR_TICKS_AND_NUMERIC_LABELS
            t.right_edge_type = AXIS_WITH_MAJOR_TICKS_ONLY
            plot_v_Star
        end
    end

    def plot_WIMP_population

        setup_data
        set_track_for_plots
        t.rescale(0.8)
        t.ylabel_scale = 0.6
        t.yaxis_numeric_label_scale = 0.6
        num_plots = 4; row = 1
        xs = d.star_Age[@track_first .. @track_last]
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 1)) do # first
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.no_xlabel
            plot_n_WIMPs
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            t.no_xlabel
            plot_cap
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            t.no_xlabel
            plot_ann
        row += 1
        end
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.top_edge_type = AXIS_HIDDEN
            plot_dNdt
        end
    end

    def plot_WIMP_distribution

        setup_data
        set_track_for_plots
        t.legend_text_dy = 1.3
        t.rescale(0.8)
        t.yaxis_numeric_label_scale = 0.6
        t.ylabel_scale = 0.7
        num_plots = 3; row = 1
        xs = d.star_Age[@track_first .. @track_last]
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 1)) do # first
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.no_xlabel
            plot_r_chi(xs)
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            t.no_xlabel
            plot_temps(xs)
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.top_edge_type = AXIS_HIDDEN
            plot_densities(xs)
        end
    end
    
end
