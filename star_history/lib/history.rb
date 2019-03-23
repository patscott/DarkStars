# history.rb

class StarHistory

    include Math
    include Tioga
    include FigureConstants
    
    attr_reader :age_Units, :age_Scale, :age_String, :profile_names, :mesh_lines_flag
    
    def initialize(dir, dir2, mesh_flag=false)
        star_history_init(dir, dir2, mesh_flag)
        return self
    end
    
    def get_num_profiles
        fname = sprintf("%s/profiles.log", @history_dir)
        file = File.open(fname)
        @profile_names = []
        while true do
            model = -1
            line = file.scanf("%s %s %s")
            if line == nil or line.length < 3
                file.close
                return @profile_names.length
            end
            age = line[0]
            mass = line[1]
            model = line[2]
            name = sprintf("model_%i.log", model)
            fullname = @profile_dir + '/' + name
            @profile_names << fullname
        end
    end
    
    def star_history_init(dir, dir2, mesh_flag=false)
        @mesh_lines_flag = mesh_flag
        @figure_maker = FigureMaker.default
        @history_dir = dir
        @profile_dir = dir2
        @track_start_param =  0.1 # fractional way from 1st to 2nd profiles
        @track_end_param = 0.0 # fractional way from last to next to last profiles
        t.def_eval_function { |str| eval(str) }

        t.save_dir = 'history_out'
        t.def_figure("H_R") { plot_H_R }
        t.def_figure("T_RHO") { plot_T_RHO }
        t.def_figure("H_R_T_RHO") { plot_H_R_T_RHO }
        t.def_figure("R_L_Ts_Tc_Rho_Psi") { plot_R_L_Ts_Tc_Rho_Psi }
        t.def_figure("Z_3a_aZ_nu_M_Core") { plot_Z_3a_aZ_nu_M_Core }
        t.def_figure("Abund_Power_Neu_Pressure") { plot_abund_power_nu_press }
        t.def_figure("Convection_and_Burning") { plot_convection_epsnuc }
        t.def_figure("Magnitudes_1") { plot_magnitudes_1 }
        t.def_figure("Magnitudes_2") { plot_magnitudes_2 }
        t.def_figure("WIMP_luminosity") {plot_WIMP_luminosity}
        t.def_figure("WIMP_population") {plot_WIMP_population}
        t.def_figure("WIMP_distribution") {plot_WIMP_distribution}
        t.def_figure("WIMP_combo") {plot_WIMP_combo}
        @num_profiles = get_num_profiles
        t.def_figure("Profile1") { plot_profile(1) } if @num_profiles >= 1
        t.def_figure("EoSProfile1") { plot_eos_profile(1) } if @num_profiles >= 1
        t.def_figure("Profile2") { plot_profile(2) } if @num_profiles >= 2
        t.def_figure("EoSProfile2") { plot_eos_profile(2) } if @num_profiles >= 2
        t.def_figure("Profile3") { plot_profile(3) } if @num_profiles >= 3
        t.def_figure("EoSProfile3") { plot_eos_profile(3) } if @num_profiles >= 3
        t.def_figure("Profile4") { plot_profile(4) } if @num_profiles >= 4
        t.def_figure("EoSProfile4") { plot_eos_profile(4) } if @num_profiles >= 4
        t.def_figure("Profile5") { plot_profile(5) } if @num_profiles >= 5
        t.def_figure("EoSProfile5") { plot_eos_profile(5) } if @num_profiles >= 5
        t.def_figure("Profile6") { plot_profile(6) } if @num_profiles >= 6
        t.def_figure("EoSProfile6") { plot_eos_profile(6) } if @num_profiles >= 6
        t.def_figure("Profile7") { plot_profile(7) } if @num_profiles >= 7
        t.def_figure("EoSProfile7") { plot_eos_profile(7) } if @num_profiles >= 7
        t.def_figure("Profile8") { plot_profile(8) } if @num_profiles >= 8
        t.def_figure("EoSProfile8") { plot_eos_profile(8) } if @num_profiles >= 8
        t.def_figure("Profile9") { plot_profile(9) } if @num_profiles >= 9
        t.def_figure("EoSProfile9") { plot_eos_profile(9) } if @num_profiles >= 9
        @profile_data = []
        @profile_names.each { |name|
            @profile_data << data = ProfileData.new
            data.read_profile(name) }
        @profiler = ProfilePlots.new(nil, nil, true)
        @have_setup_data = false
        @margin = 0.1
        
        t.def_enter_page_function { enter_page }
            
    end
    
    def enter_page
        t.page_setup(11*72/2,8.5*72/2)
        t.set_frame_sides(0.15,0.85,0.85,0.15) # left, right, top, bottom in page coords        
    end
    
    def plot_profile(num)
        data = @profile_data[num-1]
        @profiler.set_profile_data(data)
        @profiler.full_profile
    end
    
    def plot_eos_profile(num)
        data = @profile_data[num-1]
        @profiler.set_profile_data(data)
        @profiler.entropy_eos
    end
    
    def t
        @figure_maker
    end

    def d
        @log
    end

    def read_PSIs(path = nil)
        path = '' if path == nil
        path = path + '/' if path.length > 0 && path[-1..-1] != '/'
        Dvector.read(path + 'psi0.data', @psi0 = [Dvector.new, Dvector.new], 2)
        Dvector.read(path + 'psi5.data', @psi5 = [Dvector.new, Dvector.new], 2)
        Dvector.read(path + 'psi100.data', @psi100 = [Dvector.new, Dvector.new], 2)
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
    
    def setup_data
        return if @have_setup_data
        @log = HistoryLogs.new(@history_dir, @mesh_lines_flag)
        read_PSIs('star_data')
        read_ZAMS('star_data')
        adjust_ages
        find_He_gap
        @have_setup_data = true
    end
    
    def adjust_ages
        initial_age = d.star_Age[0]  #for zeroing the age after restart
        d.star_Age.sub!(initial_age) #for zeroing the age after restart
        age = d.star_Age[-1]
        if age >= 1e9
            d.star_Age.times!(1e-9)
            d.profile_Ages.times!(1e-9)
            @age_Units = "Gyr"
            @age_Scale = 1e9
            @age_xlabel = 'Age (Gyr)'
            @age_String = '$10^9$'
        elsif age >= 1e6
            d.star_Age.times!(1e-6)
            d.profile_Ages.times!(1e-6)
            @age_Units = "Myr"
            @age_Scale = 1e6
            @age_xlabel = 'Age (Myr)'
            @age_String = '$10^6$'
        elsif age >= 1e3
            d.star_Age.times!(1e-3)
            d.profile_Ages.times!(1e-3)
            @age_Units = "Kyr"
            @age_Scale = 1e3
            @age_xlabel = 'Age (Kyr)'
            @age_String = '$10^3$'
        else
            @age_Units = "yr"
            @age_Scale = 1
            @age_xlabel = 'Age beyond 500\,Myr (yr)'
            @age_String = ''
        end
      @age_format = (d.star_Age[-1] > 9.99)? '\sffamily %0.2f' : '\sffamily %0.3f'
      @mass_format = (d.star_Mass[0] > 9.99)? '\sffamily %0.2f' : '\sffamily %0.3f'
    end
    
    def find_He_gap # may have a gap in tracks if had a helium flash
        @he_burn_start = d.net_He_power.where_gt(0.001)
        if @he_burn_start == nil || @he_burn_start >= d.num_models-1
            @gap_in_tracks = false
        else
            density_jumps = d.log_center_Density[@he_burn_start .. -2] - d.log_center_Density[@he_burn_start+1 .. -1]
            jumps_at = density_jumps.where_gt(0.9) # large drop in density signals gap in tracks
            if jumps_at == nil
                @gap_in_tracks = false
            else
                @helium_gap_start = @he_burn_start + 1 + jumps_at
                @gap_in_tracks = (@helium_gap_start != nil && @helium_gap_start < d.num_models)
            end
        end
    end
    
    def set_track_for_plots
        if d.num_profiles < 4
            @track_first = 0
            @track_last = -1
        else
            start_age = (1-@track_start_param)*d.profile_Ages[0]+@track_start_param*d.profile_Ages[1]
            @track_first = d.star_Age.where_closest(start_age)
            end_age = (1-@track_end_param)*d.profile_Ages[-1]+@track_end_param*d.profile_Ages[-2]
            @track_last = d.star_Age.where_closest(end_age)
        end
        @track_last -= 1 if @track_last+1 == d.star_Age[-1] # don't use the very last model
        @start_age = d.star_Age[@track_first]
        @start_age = d.star_Age[@track_first]
        @end_age = d.star_Age[@track_last]
        age_range = @end_age - @start_age
        @end_age += 0.01 * age_range
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
    
    def mark_profiles_on_x
        base = t.bounds_bottom
        height = 0.8 * t.default_text_height_dy
        if t.yaxis_reversed
            top = base - height
        else
            top = base + height
        end
        t.line_type = Line_Type_Solid
        t.stroke_color = Green
        d.profile_Ages.each { |age| t.stroke_line(age, base, age, top) }
    end
    
    def mark_spot(x, y, color = Red)
        t.show_marker('x' => x, 'y' => y, 'marker' => Bullet, 'scale' => 0.5, 'color' => color);
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
    
    def show_box_labels(title, xlabel=nil, ylabel=nil)
        if title != nil
            t.show_title(title); t.no_title
        end
        if xlabel != nil
            t.show_xlabel(xlabel); t.no_xlabel
        end
        if ylabel != nil
            t.show_ylabel(ylabel); t.no_ylabel
        end
    end
    
    def background
        t.fill_color = OldLace
        t.fill_frame
    end
    
end

