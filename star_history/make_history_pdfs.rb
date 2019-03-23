    require 'Tioga/tioga'
    include Tioga
    require 'history.rb'

    t = FigureMaker.default
    t.num_figures.times { |i| t.make_preview_pdf(i) }
    
