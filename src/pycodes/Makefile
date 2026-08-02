# Makefile for Approximation Algorithms - LaTeX Book

BOOK = book
TEX = $(BOOK).tex
PDF = $(BOOK).pdf

# Main targets
all: $(PDF)

$(PDF): $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)
	bibtex $(BOOK) || true
	pdflatex -interaction=nonstopmode $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)

# Quick compile (no bibtex)
quick: $(TEX)
	pdflatex -interaction=nonstopmode $(TEX)

# Clean auxiliary files
clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.lof *.lot *.idx *.ind *.ilg *.synctex.gz *.fdb_latexmk *.fls

# Clean everything including PDF
distclean: clean
	rm -f $(PDF)

# Run Python demos
demo:
	python3 main.py

# Vazirani chapters
ch1:
	python3 intro.py

ch2:
	python3 set_cover.py

ch3:
	python3 steiner_tsp.py

ch4:
	python3 multiway_kcut.py

ch5:
	python3 kcenter.py

ch6:
	python3 feedback_vertex_set.py

ch7:
	python3 shortest_superstring.py

ch8:
	python3 knapsack.py

ch9:
	python3 bin_packing.py

ch10:
	python3 makespan.py

ch11:
	python3 euclidean_tsp.py

ch12:
	python3 lp_algorithms.py

ch15:
	python3 weighted_vertex_cover_pd.py

ch16:
	python3 randomized_rounding.py

ch17:
	python3 chernoff_bounds.py

ch18:
	python3 sdp_maxcut.py

ch19:
	python3 multiway_rounding.py

ch21:
	python3 steiner_forest.py

ch22:
	python3 steiner_network.py

ch23:
	python3 primal_dual_fvs.py

ch24:
	python3 facility_location.py

ch25:
	python3 kmedian.py

ch26:
	python3 sdp_sat.py

ch27:
	python3 shortest_vector.py

ch28:
	python3 counting_problems.py

ch29:
	python3 hardness_approx.py

ch30:
	python3 tree_multicut.py

ch_unrelated:
	python3 unrelated_scheduling.py

ch_multicut:
	python3 multicut_general.py

ch_sparsest:
	python3 sparsest_cut.py

ch_multicut_trees:
	python3 multicut_trees.py

# Williamson & Shmoys chapters
ws_ch2:
	python3 edge_coloring.py

ws_ch4:
	python3 prize_collecting_steiner.py

ws_ch7:
	python3 shortest_st_path_pd.py
	python3 generalized_steiner.py
	python3 min_knapsack_pd.py

ws_ch8:
	python3 balanced_cuts.py
	python3 tree_metrics.py
	python3 buy_at_bulk.py

ws_ch9:
	python3 facility_location_ls.py
	python3 kmedian_ls.py

ws_ch11:
	python3 iterated_rounding.py

ws_ch12:
	python3 steiner_tree_rr.py
	python3 rent_or_buy.py

ws_ch15:
	python3 oblivious_routing.py

ws_all:
	python3 edge_coloring.py
	python3 prize_collecting_steiner.py
	python3 shortest_st_path_pd.py
	python3 generalized_steiner.py
	python3 min_knapsack_pd.py
	python3 balanced_cuts.py
	python3 tree_metrics.py
	python3 buy_at_bulk.py
	python3 facility_location_ls.py
	python3 kmedian_ls.py
	python3 iterated_rounding.py
	python3 steiner_tree_rr.py
	python3 rent_or_buy.py
	python3 oblivious_routing.py

.PHONY: all quick clean distclean demo ch1 ch2 ch3 ch4 ch5 ch6 ch7 ch8 ch9 ch10 ch11 ch12 ch15 ch16 ch17 ch18 ch19 ch21 ch22 ch23 ch24 ch25 ch26 ch27 ch28 ch29 ch30 ws_ch2 ws_ch4 ws_ch7 ws_ch8 ws_ch9 ws_ch11 ws_ch12 ws_ch15 ws_all
