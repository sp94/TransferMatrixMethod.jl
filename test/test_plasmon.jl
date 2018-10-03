# Compare results to TMM code by Steve Byrnes, https://github.com/sbyrnes321/tmm
benchmark_source = """
%pylab
from tmm import coh_tmm
lam0 = 632 # nanometers
d_list = [np.inf, 50.0, np.inf] # nanometers
n_list = [1.5, sqrt(-16+0.5j), 1.0]
r_list, t_list = [], []
angles = np.radians(np.arange(0,41,5))
benchmark = []
for angle in angles:
    results = coh_tmm("p", n_list, d_list, angle, lam0)
    benchmark.append((results["r"],results["t"]))
"""
benchmark = [
  [0.7224288187574386+0.6571932371176172im, 0.11468484403568181-0.14623230040908508im],
  [0.7204635511861754+0.6591341476267213im, 0.11599210227582156-0.14682107962178634im],
  [0.7144800187884609+0.6649664158476095im, 0.12008751331405022-0.14862649732876415im],
  [0.7042093713302134+0.6747115862613352im, 0.12755848583031543-0.15177068619497170im],
  [0.6891862566362971+0.6883771328810835im, 0.13966136364923740-0.15647031245448442im],
  [0.6687334654105929+0.7058906756737224im, 0.15902186563221400-0.16302743602138364im],
  [0.6419975224290142+0.7269246115729855im, 0.19177695385481600-0.17159822956482856im],
  [0.6083788277506559+0.7503115166101243im, 0.25631461756401660-0.17981608390709597im],
  [0.5735314386928214+0.7712740072849034im, 0.44940400442051540-0.13417608040807666im]
]

k0 = 2pi / 632 # 1/nanometers
d_list = [0, 50, 0] # nanometers
epr_list = [2.25, -16+0.5im, 1.0]
mur_list = ones(length(epr_list))
g = Geometry(d_list, epr_list, mur_list)
results = []
for θ in (0:5:40)*pi/180
	res = tmm(g, k0, θ, 0)
	push!(results, [res.ref.p, res.trn.p])
end

@test isapprox(results, benchmark)
