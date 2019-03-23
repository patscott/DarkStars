    require 'Tioga/tioga'
    include Tioga
    require 'load.rb'

    file = File.open('which_profile_number.txt') # a number from 1 to 6 or 7 or however many profiles there are
    num = file.gets.chop
    file.close
    name = GetModelName.getname(num.to_i)
    puts name
    puts " "
    ProfilePlots.new(name)
    t = FigureMaker.default
    t.num_figures.times { |i| t.make_preview_pdf(i) }
    
