    require 'Tioga/tioga'
    include Tioga
    load 'PatPlots.rb'
    load 'PatPlots2.rb'

    t = FigureMaker.default
    PatPlots.new
    PatPlots2.new
    t.num_figures.times { |i| t.make_preview_pdf(i) }
    
