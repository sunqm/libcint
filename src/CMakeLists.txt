target_sources(Cint
  PRIVATE
	c2f.c
	cart2sph.c
	cint1e.c
	cint2e.c
	cint_bas.c
  	fblas.c
	g1e.c
	g2e.c
	misc.c
	optimizer.c
  	rys_roots.c
  	cint2c2e.c
	g2c2e.c
	cint3c2e.c
	g3c2e.c
  	cint3c1e.c
	g3c1e.c
	breit.c
	cint4c1e.c
	g4c1e.c
  	cint1e_a.c
	cint3c1e_a.c
)

if(WITH_COULOMB_ERF)
	target_sources(Cint g2e_coulerf.c cint2e_coulerf.c)

endif()

if(WITH_F12)
	target_sources(Cint g2e_f12.c stg_roots.c cint2e_f12.c)
endif()

if(WITH_GTG)
	target_sources(Cint g2e_gtg.c cint2e_gtg.c cint3c2e_gtg.c cint2c2e_gtg.c)
endif()

add_subdirectory(autocode)
