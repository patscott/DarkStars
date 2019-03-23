# plots.rb

class StarHistory

    include Math
    include FigureConstants
    
    def mark_line_at_loc(xs, ys, loc, color=Red)
        value = ys[loc]
        x = xs[loc]
        mark_spot(x, value, color)
        abs_val = value.abs
        if abs_val >= 100
            fmt = '\sffamily$\;\;$%0.1f'
        elsif abs_val >= 10
            fmt = '\sffamily$\;\;$%0.2f'
        else
            fmt = '\sffamily$\;\;$%0.3f'
        end
        t.show_label('text' => sprintf(fmt, value), 'x' => x, 'y' => value,
            'alignment' => ALIGNED_AT_MIDHEIGHT,
            'color' => Black, 'justification' => LEFT_JUSTIFIED, 'scale' => 0.9)
    end
    
    def mark_last_value(xs, ys)
        mark_spot(xs[-1], ys[-1])
        mark_line_at_loc(xs, ys, -1)
    end
    
    def plot_boundaries(xs,ys,margin,ymin=nil,ymax=nil,ymargin=nil)
        extra_on_left = 0.1
        extra_on_right = 1.5
        extra_on_top = 2.5
        extra_on_bottom = 1.5
        xmin = xs.min
        xmax = xs.max
        ymin = ys.min if ymin == nil
        ymax = ys.max if ymax == nil
        width = (xmax == xmin)? 1 : xmax - xmin
        height = (ymax == ymin)? 1 : ymax - ymin
        xmargin = margin
        ymargin = margin if ymargin == nil
        left_boundary = xmin - extra_on_left * xmargin * width
        right_boundary = xmax + extra_on_right * xmargin * width
        top_boundary = ymax + extra_on_top * ymargin * height
        bottom_boundary = ymin - extra_on_bottom * ymargin * height
        return [ left_boundary, right_boundary, top_boundary, bottom_boundary ]
    end
    
    def plot_one_value(ys, ylabel, ymin_limit = nil)
        xs = d.star_Age[@track_first .. @track_last]
        ys = ys[@track_first .. @track_last]
        t.ylabel_scale = 0.8
        t.show_ylabel(ylabel); t.no_ylabel
        t.yaxis_numeric_label_scale = 0.6
        ymin = ys.min
        ymin = ymin_limit if ymin_limit != nil && ymin < ymin_limit
        t.show_plot('boundaries' => plot_boundaries(xs,ys,@margin,ymin)) do
            background
            mark_profiles_on_x
            t.stroke_color = Black
            stroke_track(xs,ys)
            mark_last_value(xs, ys)
        end
    end
    
    def plot_R
        plot_one_value(d.log_Radius, '$\log$ $\mathrm{R/R_{\odot}}$')
    end
    
    def plot_L
        plot_one_value(d.log_Luminosity, '$\log$ $\mathrm{L/L_{\odot}}$')
    end
    
    def plot_Ts
        plot_one_value(d.log_surface_Temp, '$\log$ T surf')
    end
    
    def plot_Tc
        plot_one_value(d.log_center_Temp, '$\log$ T cntr')
    end
    
    def plot_Rho
        plot_one_value(d.log_center_Density, '$\log$ $\rho$ cntr')
    end
    
    def plot_Psi
        plot_one_value(d.center_Degeneracy, '$\Psi$ cntr')
    end
    
    def plot_Z
        plot_one_value(1.0 - d.center_H - d.center_He, 'Z cntr')
    end
    
    def plot_3alpha
        plot_one_value(d.power_3_alpha.safe_log10, '$\log$ $3 \alpha$ $\mathrm{L/L_{\odot}}$', -6)
    end
    
    def plot_alpha_Metal
        plot_one_value((d. power_He_burn - d.power_3_alpha).safe_log10, '$\log$ $\alpha$+metal', -6)
    end
    
    def plot_Neutrino_Loss
        plot_one_value(d.power_Neutrinos.safe_log10, '$\log$ $\nu$ loss', -6)
    end
    
    def plot_Mass
        plot_one_value(d.star_Mass, '$\mathrm{M/M_{\odot}}$')
    end
    
    def plot_Core
        plot_one_value(d.mass_He_Core, 'Core $\mathrm{M/M_{\odot}}$')
    end
    
    def plot_R_L_Ts_Tc_Rho_Psi
        setup_data
        set_track_for_plots
        t.rescale(0.6)
        num_plots = 6; row = 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 1)) do # first
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            plot_R
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_L
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_Ts
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_Tc
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_Rho
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => num_plots)) do # last
            t.top_edge_type = AXIS_HIDDEN
            plot_Psi
            t.show_xlabel(@age_xlabel)
        end
    end
    
    def plot_Z_3a_aZ_nu_M_Core
        setup_data
        set_track_for_plots
        t.rescale(0.6)
        num_plots = 6; row = 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 1)) do # first
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            plot_Z
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_3alpha
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_alpha_Metal
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_Neutrino_Loss
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_Mass
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => num_plots)) do # last
            t.top_edge_type = AXIS_HIDDEN
            plot_Core
            t.show_xlabel(@age_xlabel)
        end
    end
    
    def plot_u_minus_b
        plot_one_value(d.umb, 'U-B')
    end
    
    def plot_b_minus_v
        plot_one_value(d.bmv, 'B-V')
    end
    
    def plot_v_minus_r
        plot_one_value(d.vmr, 'V-R')
    end
    
    def plot_v_minus_i
        plot_one_value(d.vmi, 'V-I')
    end
    
    def plot_v_minus_k
        plot_one_value(d.vmk, 'V-K')
    end
    
    def plot_r_minus_i
        plot_one_value(d.rmi, 'R-I')
    end
    
    def plot_i_minus_k
        plot_one_value(d.imk, 'I-K')
    end
    

    def plot_j_minus_h
        plot_one_value(d.jmh, 'J-H')
    end
    
    def plot_h_minus_k
        plot_one_value(d.hmk, 'H-K')
    end
    
    def plot_k_minus_l
        plot_one_value(d.kml, 'K-L')
    end
    
    def plot_j_minus_k
        plot_one_value(d.jmk, 'J-K')
    end
    
    def plot_j_minus_l
        plot_one_value(d.jml, 'J-L')
    end
    
    def plot_j_minus_lp
        plot_one_value(d.jmlp, 'J-LP')
    end
    
    def plot_k_minus_m
        plot_one_value(d.kmm, 'K-M')
    end
    
    def plot_v_mag
        plot_one_value(d.bol - d.bcv, 'V mag')
    end
    

    def plot_magnitudes_1
        setup_data
        set_track_for_plots
        t.rescale(0.6)
        num_plots = 7; row = 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 1)) do # first
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            plot_v_mag
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_u_minus_b
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_b_minus_v
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_v_minus_r
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_v_minus_i
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_v_minus_k
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => num_plots)) do # last
            t.top_edge_type = AXIS_HIDDEN
            plot_r_minus_i
            t.show_xlabel(@age_xlabel)
        end
    end
    
    
    def plot_magnitudes_2
        setup_data
        set_track_for_plots
        t.rescale(0.6)
        num_plots = 7; row = 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => 1)) do # first
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            plot_i_minus_k
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_j_minus_h
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_h_minus_k
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_k_minus_l
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_j_minus_k
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => row)) do 
            t.xaxis_type = AXIS_WITH_TICKS_ONLY
            t.top_edge_type = AXIS_HIDDEN
            plot_j_minus_l
        end
        row += 1
        t.subplot(t.row_margins('num_rows' => num_plots, 'row' => num_plots)) do # last
            t.top_edge_type = AXIS_HIDDEN
            plot_k_minus_m
            t.show_xlabel(@age_xlabel)
        end
    end
    
end
