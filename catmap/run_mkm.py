from catmap import ReactionModel

mkm_file = 'furfural_hydro_vacuum.mkm'
model = ReactionModel(setup_file=mkm_file)
model.output_variables += ['rate', 'directional_rates', 'production_rate','coverage']
model.run()

from catmap import analyze
vm = analyze.VectorMap(model)


vm.plot_variable = 'rate'
vm.log_scale = True
vm.min = 1e-10
vm.max = 1e6
fig2 = vm.plot(save='rate.png')

vm.plot_variable = 'coverage'
vm.min = 0
vm.max = 1
vm.log_scale = False
fig3 = vm.plot(save='coverage.png')

sa = analyze.ScalingAnalysis(model)
sa.plot(save='scaling.pdf')

vm = analyze.VectorMap(model)
vm.plot_variable = 'production_rate' #tell the model which output to plot
vm.log_scale = True #rates should be plotted on a log-scale
vm.min = 1e-20 #minimum rate to plot
vm.max = 1e8 #maximum rate to plot
vm.threshold = 1e-25 #anything below this is considered to be 0
vm.subplots_adjust_kwargs = {'left':0.2,'right':0.8,'bottom':0.15}
vm.plot(save='production_rate.png')

