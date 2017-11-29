# DISCRETA makefile
# Anton Betten 1995
#
# 

#CONF=CONF_LONDON97
#CFLAGS=-O2 -D$(CONF) -DSYSTEM_UNIX -Wunused
# -Wunused
#CC=g++

A2PSOPTS=-nP
LPS=rsh btm2xd lp -dps
#LPS= > /dev/null

INCLUDES = -I..

OBJPATH = ../obj/$(DISCRETA_ARCH)

OBJ = $(OBJPATH)/DISCRETA

MODULES = \
		$(OBJ)/discreta.o \
		$(OBJ)/in.o \
		$(OBJ)/in1.o \
		$(OBJ)/io.o \
		$(OBJ)/io0.o \
		$(OBJ)/io1.o \
		$(OBJ)/io2.o \
		$(OBJ)/iof.o \
		$(OBJ)/nu.o \
		$(OBJ)/os.o \
		$(OBJ)/bch.o \
		$(OBJ)/mindist.o \
		$(OBJ)/divs.o \
		$(OBJ)/list.o \
		$(OBJ)/ma.o \
		$(OBJ)/ma_geo.o \
		$(OBJ)/unip.o \
		$(OBJ)/vec.o  \
	$(OBJ)/bt.o \
	$(OBJ)/bt_key.o \
	$(OBJ)/db.o \
	$(OBJ)/db1.o \
	$(OBJ)/dp1.o \
	$(OBJ)/dp_html.o \
	$(OBJ)/generators.o \
	$(OBJ)/conway.o \
	$(OBJ)/dc_draw.o \
	$(OBJ)/dp.o \
	$(OBJ)/graphical.o \
	$(OBJ)/higman_sims.o \
	$(OBJ)/intersection.o \
	$(OBJ)/intersection_aijk.o \
	$(OBJ)/kramer_mesner.o \
	$(OBJ)/ladder.o \
	$(OBJ)/ladder2.o \
	$(OBJ)/ladder_info.o \
	$(OBJ)/mathieu.o \
	$(OBJ)/parameter.o \
	$(OBJ)/plesken.o \
	$(OBJ)/report.o \
		$(OBJ)/cp.o \
		$(OBJ)/fga.o \
		$(OBJ)/fga_gen.o \
		$(OBJ)/kontext.o \
		$(OBJ)/lb.o \
		$(OBJ)/lb_base_change.o \
		$(OBJ)/mon.o \
		$(OBJ)/part.o \
		$(OBJ)/perm.o \
		$(OBJ)/perm2.o \
		$(OBJ)/perm_grp.o \
		$(OBJ)/perm_grp2.o \
		$(OBJ)/poly.o \
		$(OBJ)/single_cosets.o \
		$(OBJ)/dcc_orderly.o \
		$(OBJ)/geo_canon.o \
		$(OBJ)/geo_data.o  \
		$(OBJ)/geo_store.o  \
		$(OBJ)/iso.o \
		$(OBJ)/iso_sub.o \
		$(OBJ)/ntdo.o \
		$(OBJ)/ntdo2.o \
		$(OBJ)/ntdo_dd.o \
		$(OBJ)/tda.o \
		$(OBJ)/tdo.o \
		$(OBJ)/bruch.o \
		$(OBJ)/gfq_nb.o \
		$(OBJ)/gfq_psu.o \
		$(OBJ)/gfq_sz.o \
		$(OBJ)/gfq_zech.o \
		$(OBJ)/lo.o \
		$(OBJ)/perm_rep.o \
		$(OBJ)/singer.o \
	$(OBJ)/epic.o \
	$(OBJ)/ged.o \
	$(OBJ)/mp.o \
	$(OBJ)/ps.o \
	$(OBJ)/tree.o \
	$(OBJ)/vbp.o \
	$(OBJ)/vdi.o \
	$(OBJ)/dimino.o \
	$(OBJ)/sgl.o \
	$(OBJ)/sgld.o \
	$(OBJ)/sglo.o \
	$(OBJ)/sgls.o \
	$(OBJ)/sgls2.o \
	$(OBJ)/ruckdeschel.o \
		$(OBJ)/aut.o \
		$(OBJ)/cf.o \
		$(OBJ)/cl_rep.o \
		$(OBJ)/cl_rep2.o \
		$(OBJ)/gt_canon.o \
		$(OBJ)/gt_col_util.o \
		$(OBJ)/gt_color.o \
		$(OBJ)/fg.o \
		$(OBJ)/fg_color.o \
		$(OBJ)/fg_direct.o \
		$(OBJ)/fg_table.o \
		$(OBJ)/fg_iso.o \
		$(OBJ)/fg_syl.o \
		$(OBJ)/fg_ext.o 

all : $(MODULES) lib;

lib: $(MODULES) ;
	ar ruvs $(OBJPATH)/libDISCRETA.a $(OBJ)/*.o

$(OBJ)/discreta.o: discreta.C ; $(CC) -c $(CFLAGS) $(INCLUDES) discreta.C -o $(OBJ)/discreta.o
$(OBJ)/in.o: in.C ; $(CC) -c $(CFLAGS) $(INCLUDES) in.C -o $(OBJ)/in.o
$(OBJ)/in1.o: in1.C ; $(CC) -c $(CFLAGS) $(INCLUDES) in1.C -o $(OBJ)/in1.o
$(OBJ)/io.o: io.C ; $(CC) -c $(CFLAGS) $(INCLUDES) io.C -o $(OBJ)/io.o
$(OBJ)/io0.o: io0.C ; $(CC) -c $(CFLAGS) $(INCLUDES) io0.C -o $(OBJ)/io0.o
$(OBJ)/io1.o: io1.C ; $(CC) -c $(CFLAGS) $(INCLUDES) io1.C -o $(OBJ)/io1.o
$(OBJ)/io2.o: io2.C ; $(CC) -c $(CFLAGS) $(INCLUDES) io2.C -o $(OBJ)/io2.o
$(OBJ)/iof.o: iof.C ; $(CC) -c $(CFLAGS) $(INCLUDES) iof.C -o $(OBJ)/iof.o
$(OBJ)/nu.o: nu.C ; $(CC) -c $(CFLAGS) $(INCLUDES) nu.C -o $(OBJ)/nu.o
$(OBJ)/os.o: os.C ; $(CC) -c $(CFLAGS) $(INCLUDES) os.C -o $(OBJ)/os.o

$(OBJ)/bch.o: bch.C ; $(CC) -c $(CFLAGS) $(INCLUDES) bch.C -o $(OBJ)/bch.o
$(OBJ)/mindist.o: mindist.C ; $(CC) -c $(CFLAGS) $(INCLUDES) mindist.C -o $(OBJ)/mindist.o

$(OBJ)/divs.o: divs.C ; $(CC) -c $(CFLAGS) $(INCLUDES) divs.C -o $(OBJ)/divs.o
$(OBJ)/list.o: list.C ; $(CC) -c $(CFLAGS) $(INCLUDES) list.C -o $(OBJ)/list.o
$(OBJ)/ma.o: ma.C ; $(CC) -c $(CFLAGS) $(INCLUDES) ma.C -o $(OBJ)/ma.o
$(OBJ)/ma_geo.o: ma_geo.C ; $(CC) -c $(CFLAGS) $(INCLUDES) ma_geo.C -o $(OBJ)/ma_geo.o
$(OBJ)/unip.o: unip.C ; $(CC) -c $(CFLAGS) $(INCLUDES) unip.C -o $(OBJ)/unip.o
$(OBJ)/vec.o: vec.C ; $(CC) -c $(CFLAGS) $(INCLUDES) vec.C -o $(OBJ)/vec.o

$(OBJ)/bt.o: bt.C ; $(CC) -c $(CFLAGS) $(INCLUDES) bt.C -o $(OBJ)/bt.o
$(OBJ)/bt_key.o: bt_key.C ; $(CC) -c $(CFLAGS) $(INCLUDES) bt_key.C -o $(OBJ)/bt_key.o
$(OBJ)/db.o: db.C ; $(CC) -c $(CFLAGS) $(INCLUDES) db.C -o $(OBJ)/db.o
$(OBJ)/db1.o: db1.C ; $(CC) -c $(CFLAGS) $(INCLUDES) db1.C -o $(OBJ)/db1.o
$(OBJ)/dp1.o: dp1.C ; $(CC) -c $(CFLAGS) $(INCLUDES) dp1.C -o $(OBJ)/dp1.o
$(OBJ)/dp_html.o: dp_html.C ; $(CC) -c $(CFLAGS) $(INCLUDES) dp_html.C -o $(OBJ)/dp_html.o
$(OBJ)/generators.o: generators.C ; $(CC) -c $(CFLAGS) $(INCLUDES) generators.C -o $(OBJ)/generators.o

$(OBJ)/conway.o: conway.C ; $(CC) -c $(CFLAGS) $(INCLUDES) conway.C -o $(OBJ)/conway.o
$(OBJ)/dc_draw.o: dc_draw.C ; $(CC) -c $(CFLAGS) $(INCLUDES) dc_draw.C -o $(OBJ)/dc_draw.o
$(OBJ)/ladder.o: ladder.C ; $(CC) -c $(CFLAGS) $(INCLUDES) ladder.C -o $(OBJ)/ladder.o
$(OBJ)/ladder2.o: ladder2.C ; $(CC) -c $(CFLAGS) $(INCLUDES) ladder2.C -o $(OBJ)/ladder2.o
$(OBJ)/ladder_info.o: ladder_info.C ; $(CC) -c $(CFLAGS) $(INCLUDES) ladder_info.C -o $(OBJ)/ladder_info.o
$(OBJ)/dp.o: dp.C ; $(CC) -c $(CFLAGS) $(INCLUDES) dp.C -o $(OBJ)/dp.o
$(OBJ)/graphical.o: graphical.C ; $(CC) -c $(CFLAGS) $(INCLUDES) graphical.C -o $(OBJ)/graphical.o
$(OBJ)/higman_sims.o: higman_sims.C ; $(CC) -c $(CFLAGS) $(INCLUDES) higman_sims.C -o $(OBJ)/higman_sims.o
$(OBJ)/intersection.o: intersection.C ; $(CC) -c $(CFLAGS) $(INCLUDES) intersection.C -o $(OBJ)/intersection.o
$(OBJ)/intersection_aijk.o: intersection_aijk.C ; $(CC) -c $(CFLAGS) $(INCLUDES) intersection_aijk.C -o $(OBJ)/intersection_aijk.o
$(OBJ)/kramer_mesner.o: kramer_mesner.C ; $(CC) -c $(CFLAGS) $(INCLUDES) kramer_mesner.C -o $(OBJ)/kramer_mesner.o
$(OBJ)/mathieu.o: mathieu.C ; $(CC) -c $(CFLAGS) $(INCLUDES) mathieu.C -o $(OBJ)/mathieu.o
$(OBJ)/parameter.o: parameter.C ; $(CC) -c $(CFLAGS) $(INCLUDES) parameter.C -o $(OBJ)/parameter.o
$(OBJ)/plesken.o: plesken.C ; $(CC) -c $(CFLAGS) $(INCLUDES) plesken.C -o $(OBJ)/plesken.o
$(OBJ)/report.o: report.C ; $(CC) -c $(CFLAGS) $(INCLUDES) report.C -o $(OBJ)/report.o

$(OBJ)/cp.o: cp.C ; $(CC) -c $(CFLAGS) $(INCLUDES) cp.C -o $(OBJ)/cp.o
$(OBJ)/fga.o: fga.C ; $(CC) -c $(CFLAGS) $(INCLUDES) fga.C -o $(OBJ)/fga.o
$(OBJ)/fga_gen.o: fga_gen.C ; $(CC) -c $(CFLAGS) $(INCLUDES) fga_gen.C -o $(OBJ)/fga_gen.o
$(OBJ)/kontext.o: kontext.C ; $(CC) -c $(CFLAGS) $(INCLUDES) kontext.C -o $(OBJ)/kontext.o
$(OBJ)/lb.o: lb.C ; $(CC) -c $(CFLAGS)  $(INCLUDES) lb.C -o $(OBJ)/lb.o
$(OBJ)/lb_base_change.o: lb_base_change.C ; $(CC) -c $(CFLAGS)  $(INCLUDES) lb_base_change.C -o $(OBJ)/lb_base_change.o
$(OBJ)/mon.o: mon.C ; $(CC) -c $(CFLAGS) $(INCLUDES) mon.C -o $(OBJ)/mon.o
$(OBJ)/part.o: part.C ; $(CC) -c $(CFLAGS) $(INCLUDES) part.C -o $(OBJ)/part.o
$(OBJ)/perm.o: perm.C ; $(CC) -c $(CFLAGS) $(INCLUDES) perm.C -o $(OBJ)/perm.o
$(OBJ)/perm2.o: perm2.C ; $(CC) -c $(CFLAGS) $(INCLUDES) perm2.C -o $(OBJ)/perm2.o
$(OBJ)/perm_grp.o: perm_grp.C ; $(CC) -c $(CFLAGS) $(INCLUDES) perm_grp.C -o $(OBJ)/perm_grp.o
$(OBJ)/perm_grp2.o: perm_grp2.C ; $(CC) -c $(CFLAGS) $(INCLUDES) perm_grp2.C -o $(OBJ)/perm_grp2.o
$(OBJ)/poly.o: poly.C ; $(CC) -c $(CFLAGS) $(INCLUDES) poly.C -o $(OBJ)/poly.o
$(OBJ)/single_cosets.o: single_cosets.C ; $(CC) -c $(CFLAGS) $(INCLUDES) single_cosets.C -o $(OBJ)/single_cosets.o

$(OBJ)/dcc_orderly.o: dcc_orderly.C ; $(CC) -c $(CFLAGS) $(INCLUDES) dcc_orderly.C -o $(OBJ)/dcc_orderly.o
$(OBJ)/geo_canon.o: geo_canon.C ; $(CC) -c $(CFLAGS) $(INCLUDES) geo_canon.C -o $(OBJ)/geo_canon.o
$(OBJ)/geo_data.o: geo_data.C ; $(CC) -c $(CFLAGS) $(INCLUDES) geo_data.C -o $(OBJ)/geo_data.o
$(OBJ)/geo_store.o: geo_store.C ; $(CC) -c $(CFLAGS) $(INCLUDES) geo_store.C -o $(OBJ)/geo_store.o
$(OBJ)/iso.o: iso.C ; $(CC) -c $(CFLAGS) $(INCLUDES) iso.C -o $(OBJ)/iso.o
$(OBJ)/iso_sub.o: iso_sub.C ; $(CC) -c $(CFLAGS) $(INCLUDES) iso_sub.C -o $(OBJ)/iso_sub.o
$(OBJ)/ntdo.o: ntdo.C ; $(CC) -c $(CFLAGS) $(INCLUDES) ntdo.C -o $(OBJ)/ntdo.o
$(OBJ)/ntdo2.o: ntdo2.C ; $(CC) -c $(CFLAGS) $(INCLUDES) ntdo2.C -o $(OBJ)/ntdo2.o
$(OBJ)/ntdo_dd.o: ntdo_dd.C ; $(CC) -c $(CFLAGS) $(INCLUDES) ntdo_dd.C -o $(OBJ)/ntdo_dd.o
$(OBJ)/tda.o: tda.C ; $(CC) -c $(CFLAGS) $(INCLUDES) tda.C -o $(OBJ)/tda.o
$(OBJ)/tdo.o: tdo.C ; $(CC) -c $(CFLAGS) $(INCLUDES) tdo.C -o $(OBJ)/tdo.o

$(OBJ)/bruch.o: bruch.C ; $(CC) -c $(CFLAGS) $(INCLUDES) bruch.C -o $(OBJ)/bruch.o
$(OBJ)/gfq_nb.o: gfq_nb.C ; $(CC) -c $(CFLAGS) $(INCLUDES) gfq_nb.C -o $(OBJ)/gfq_nb.o
$(OBJ)/gfq_psu.o: gfq_psu.C ; $(CC) -c $(CFLAGS) $(INCLUDES) gfq_psu.C -o $(OBJ)/gfq_psu.o
$(OBJ)/gfq_sz.o: gfq_sz.C ; $(CC) -c $(CFLAGS) $(INCLUDES) gfq_sz.C -o $(OBJ)/gfq_sz.o
$(OBJ)/gfq_zech.o: gfq_zech.C ; $(CC) -c $(CFLAGS) $(INCLUDES) gfq_zech.C -o $(OBJ)/gfq_zech.o
$(OBJ)/lo.o: lo.C ; $(CC) -c $(CFLAGS) $(INCLUDES) lo.C -o $(OBJ)/lo.o
$(OBJ)/perm_rep.o: perm_rep.C ; $(CC) -c $(CFLAGS) $(INCLUDES) perm_rep.C -o $(OBJ)/perm_rep.o
$(OBJ)/singer.o: singer.C ; $(CC) -c $(CFLAGS) $(INCLUDES) singer.C -o $(OBJ)/singer.o

$(OBJ)/epic.o: epic.C ; $(CC) -c $(CFLAGS) $(INCLUDES) epic.C -o $(OBJ)/epic.o
$(OBJ)/ged.o: ged.C ; $(CC) -c $(CFLAGS) $(INCLUDES) ged.C -o $(OBJ)/ged.o
$(OBJ)/mp.o: mp.C ; $(CC) -c $(CFLAGS) $(INCLUDES) mp.C -o $(OBJ)/mp.o
$(OBJ)/ps.o: ps.C ; $(CC) -c $(CFLAGS) $(INCLUDES) ps.C -o $(OBJ)/ps.o
$(OBJ)/tree.o: tree.C ; $(CC) -c $(CFLAGS) $(INCLUDES) tree.C -o $(OBJ)/tree.o
$(OBJ)/vbp.o: vbp.C ; $(CC) -c $(CFLAGS) $(INCLUDES) vbp.C -o $(OBJ)/vbp.o
$(OBJ)/vdi.o: vdi.C ; $(CC) -c $(CFLAGS) $(INCLUDES) vdi.C -o $(OBJ)/vdi.o

$(OBJ)/dimino.o: dimino.C ; $(CC) -c $(CFLAGS) $(INCLUDES) dimino.C -o $(OBJ)/dimino.o
$(OBJ)/sgl.o: sgl.C ; $(CC) -c $(CFLAGS) $(INCLUDES) sgl.C -o $(OBJ)/sgl.o
$(OBJ)/sgld.o: sgld.C ; $(CC) -c $(CFLAGS) $(INCLUDES) sgld.C -o $(OBJ)/sgld.o
$(OBJ)/sglo.o: sglo.C ; $(CC) -c $(CFLAGS) $(INCLUDES) sglo.C -o $(OBJ)/sglo.o
$(OBJ)/sgls.o: sgls.C ; $(CC) -c $(CFLAGS) $(INCLUDES) sgls.C -o $(OBJ)/sgls.o
$(OBJ)/sgls2.o: sgls2.C ; $(CC) -c $(CFLAGS) $(INCLUDES) sgls2.C -o $(OBJ)/sgls2.o

$(OBJ)/ruckdeschel.o: ruckdeschel.C ; $(CC) -c $(CFLAGS) $(INCLUDES) ruckdeschel.C -o $(OBJ)/ruckdeschel.o

$(OBJ)/aut.o: aut.C aut_dimino.C aut_grid.C aut_init.C aut_util.C ; $(CC) -c $(CFLAGS) $(INCLUDES) aut.C -o $(OBJ)/aut.o
$(OBJ)/cf.o: cf.C ; $(CC) -c $(CFLAGS) $(INCLUDES) cf.C -o $(OBJ)/cf.o
$(OBJ)/cl_rep.o: cl_rep.C ; $(CC) -c $(CFLAGS) $(INCLUDES) cl_rep.C -o $(OBJ)/cl_rep.o
$(OBJ)/cl_rep2.o: cl_rep2.C ; $(CC) -c $(CFLAGS) $(INCLUDES) cl_rep2.C -o $(OBJ)/cl_rep2.o
$(OBJ)/gt_canon.o: gt_canon.C gt_canon_dimino.C gt_canon_grid.C ; $(CC) -c $(CFLAGS) $(INCLUDES) gt_canon.C -o $(OBJ)/gt_canon.o
$(OBJ)/gt_col_util.o: gt_col_util.C ; $(CC) -c $(CFLAGS) $(INCLUDES) gt_col_util.C -o $(OBJ)/gt_col_util.o
$(OBJ)/gt_color.o: gt_color.C ; $(CC) -c $(CFLAGS) $(INCLUDES) gt_color.C -o $(OBJ)/gt_color.o
$(OBJ)/fg.o: fg.C ; $(CC) -c $(CFLAGS) $(INCLUDES) fg.C -o $(OBJ)/fg.o
$(OBJ)/fg_color.o: fg_color.C ; $(CC) -c $(CFLAGS) $(INCLUDES) fg_color.C -o $(OBJ)/fg_color.o
$(OBJ)/fg_direct.o: fg_direct.C ; $(CC) -c $(CFLAGS) $(INCLUDES) fg_direct.C -o $(OBJ)/fg_direct.o
$(OBJ)/fg_table.o: fg_table.C ; $(CC) -c $(CFLAGS) $(INCLUDES) fg_table.C -o $(OBJ)/fg_table.o
$(OBJ)/fg_iso.o: fg_iso.C ; $(CC) -c $(CFLAGS) $(INCLUDES) fg_iso.C -o $(OBJ)/fg_iso.o
$(OBJ)/fg_syl.o: fg_syl.C ; $(CC) -c $(CFLAGS) $(INCLUDES) fg_syl.C -o $(OBJ)/fg_syl.o
$(OBJ)/fg_ext.o: fg_ext.C ; $(CC) -c $(CFLAGS) $(INCLUDES) fg_ext.C -o $(OBJ)/fg_ext.o

ps1: ;
	a2ps $(A2PSOPTS) discreta.h | $(LPS)
	a2ps $(A2PSOPTS) discreta_types.h | $(LPS)
	a2ps $(A2PSOPTS) discreta.C | $(LPS)
	a2ps $(A2PSOPTS) in.h | $(LPS)
	a2ps $(A2PSOPTS) in.C | $(LPS)
	a2ps $(A2PSOPTS) in1.C | $(LPS)
	a2ps $(A2PSOPTS) io.C | $(LPS)
	a2ps $(A2PSOPTS) io0.C | $(LPS)
	a2ps $(A2PSOPTS) io1.C | $(LPS)
	a2ps $(A2PSOPTS) io2.C | $(LPS)
	a2ps $(A2PSOPTS) iof.h | $(LPS)
	a2ps $(A2PSOPTS) iof.C | $(LPS)
	a2ps $(A2PSOPTS) nu.C | $(LPS)
	a2ps $(A2PSOPTS) os.C | $(LPS)
	a2ps $(A2PSOPTS) divs.h | $(LPS)
	a2ps $(A2PSOPTS) divs.C | $(LPS)
	a2ps $(A2PSOPTS) list.C | $(LPS)
	a2ps $(A2PSOPTS) ma.h | $(LPS)
	a2ps $(A2PSOPTS) ma.C | $(LPS)
	a2ps $(A2PSOPTS) ma_geo.C | $(LPS)
	a2ps $(A2PSOPTS) vec.h | $(LPS)
	a2ps $(A2PSOPTS) vec.C | $(LPS)
	a2ps $(A2PSOPTS) unip.h | $(LPS)
	a2ps $(A2PSOPTS) unip.C | $(LPS)
	a2ps $(A2PSOPTS) bch.h | $(LPS)
	a2ps $(A2PSOPTS) bch.C | $(LPS)
	a2ps $(A2PSOPTS) mindist.C | $(LPS)

ps2: ;
	a2ps $(A2PSOPTS) db.h | $(LPS)
	a2ps $(A2PSOPTS) bt.C | $(LPS)
	a2ps $(A2PSOPTS) bt_key.C | $(LPS)
	a2ps $(A2PSOPTS) db.C | $(LPS)
	a2ps $(A2PSOPTS) db1.C | $(LPS)
	a2ps $(A2PSOPTS) dp.h | $(LPS)
	a2ps $(A2PSOPTS) dp.C | $(LPS)
	a2ps $(A2PSOPTS) dp1.h | $(LPS)
	a2ps $(A2PSOPTS) dp1.C | $(LPS)
	a2ps $(A2PSOPTS) dp_html.C | $(LPS)

ps3: ;
	a2ps $(A2PSOPTS) ladder.h | $(LPS)
	a2ps $(A2PSOPTS) kramer_mesner.C | $(LPS)
	a2ps $(A2PSOPTS) li.C | $(LPS)
	a2ps $(A2PSOPTS) ladder.C | $(LPS)
	a2ps $(A2PSOPTS) ladder2.C | $(LPS)
	a2ps $(A2PSOPTS) graphical.C | $(LPS)
	a2ps $(A2PSOPTS) plesken.C | $(LPS)
	a2ps $(A2PSOPTS) intersection.C | $(LPS)
	a2ps $(A2PSOPTS) intersection_aijk.C | $(LPS)
	a2ps $(A2PSOPTS) parameter.C | $(LPS)
	a2ps $(A2PSOPTS) report.C | $(LPS)
	a2ps $(A2PSOPTS) mathieu.C | $(LPS)
	a2ps $(A2PSOPTS) conway.C | $(LPS)
	a2ps $(A2PSOPTS) higman_sims.C | $(LPS)
	a2ps $(A2PSOPTS) dc_draw.h | $(LPS)
	a2ps $(A2PSOPTS) dc_draw.C | $(LPS)

ps4: ;
	a2ps $(A2PSOPTS) cp.h | $(LPS)
	a2ps $(A2PSOPTS) cp.C | $(LPS)
	a2ps $(A2PSOPTS) fga.h | $(LPS)
	a2ps $(A2PSOPTS) fga.C | $(LPS)
	a2ps $(A2PSOPTS) fga_gen.h | $(LPS)
	a2ps $(A2PSOPTS) fga_gen.C | $(LPS)
	a2ps $(A2PSOPTS) lb.h | $(LPS)
	a2ps $(A2PSOPTS) lb.C | $(LPS)
	a2ps $(A2PSOPTS) lb_base_change.C | $(LPS)
	a2ps $(A2PSOPTS) mon.C | $(LPS)
	a2ps $(A2PSOPTS) part.h | $(LPS)
	a2ps $(A2PSOPTS) part.C | $(LPS)
	a2ps $(A2PSOPTS) perm.h | $(LPS)
	a2ps $(A2PSOPTS) perm.C | $(LPS)
	a2ps $(A2PSOPTS) perm2.C | $(LPS)
	a2ps $(A2PSOPTS) perm_grp.C | $(LPS)
	a2ps $(A2PSOPTS) perm_grp2.C | $(LPS)
	a2ps $(A2PSOPTS) poly.h | $(LPS)
	a2ps $(A2PSOPTS) poly.C | $(LPS)
	a2ps $(A2PSOPTS) single_cosets.C | $(LPS)

ps5: ;
	a2ps $(A2PSOPTS) bruch.C | $(LPS)
	a2ps $(A2PSOPTS) gfq.h | $(LPS)
	a2ps $(A2PSOPTS) gfq_nb.C | $(LPS)
	a2ps $(A2PSOPTS) gfq_psu.C | $(LPS)
	a2ps $(A2PSOPTS) gfq_sz.C | $(LPS)
	a2ps $(A2PSOPTS) gfq_zech.C | $(LPS)
	a2ps $(A2PSOPTS) perm_rep.C | $(LPS)
	a2ps $(A2PSOPTS) singer.h | $(LPS)
	a2ps $(A2PSOPTS) singer.C | $(LPS)
	a2ps $(A2PSOPTS) vdi.C | $(LPS)
	a2ps $(A2PSOPTS) epic.C | $(LPS)
	a2ps $(A2PSOPTS) mp.C | $(LPS)
	a2ps $(A2PSOPTS) ps.C | $(LPS)
	a2ps $(A2PSOPTS) ged.h | $(LPS)
	a2ps $(A2PSOPTS) ged.C | $(LPS)
	a2ps $(A2PSOPTS) tree.h | $(LPS)
	a2ps $(A2PSOPTS) tree.C | $(LPS)
	a2ps $(A2PSOPTS) dimino.C | $(LPS)
	a2ps $(A2PSOPTS) sgl.C | $(LPS)
	a2ps $(A2PSOPTS) sgld.C | $(LPS)
	a2ps $(A2PSOPTS) sglo.C | $(LPS)
	a2ps $(A2PSOPTS) sgls.C | $(LPS)
	a2ps $(A2PSOPTS) sgls2.C | $(LPS)

#	a2ps $(A2PSOPTS) ruckdeschel.C | $(LPS)



ps6: ;
	a2ps $(A2PSOPTS) geo.h | $(LPS)
	a2ps $(A2PSOPTS) dcc_orderly.C | $(LPS)
	a2ps $(A2PSOPTS) geo_canon.h | $(LPS)
	a2ps $(A2PSOPTS) geo_canonP.h | $(LPS)
	a2ps $(A2PSOPTS) geo_canon.C | $(LPS)
	a2ps $(A2PSOPTS) geo_data.h | $(LPS)
	a2ps $(A2PSOPTS) geo_data.C | $(LPS)
	a2ps $(A2PSOPTS) geo_store.h | $(LPS)
	a2ps $(A2PSOPTS) geo_store.C | $(LPS)
	a2ps $(A2PSOPTS) tdo.C | $(LPS)
	a2ps $(A2PSOPTS) ntdo.C | $(LPS)
	a2ps $(A2PSOPTS) ntdo2.C | $(LPS)
	a2ps $(A2PSOPTS) ntdo_dd.C | $(LPS)
	a2ps $(A2PSOPTS) tda.C | $(LPS)



clean: ; - rm *.o *.a $(OBJ)/*.o  $(OBJ)/*.a 

tex: ;
	grepdocu -no_underscore_translation -outprefix ../obj/TEX/ *.C *.h









