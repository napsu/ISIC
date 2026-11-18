#*************************************************************************
#*                                                                       *
#*     Makefile for BigClust, the stochastic nonsmooth clustering        *
#*     software by Napsu Karmitsa 2024.                                  *
#*                                                                       *
#*     The work was financially supported by the Research Council of     *
#*     Finland (Project No. #345804 and #345805).                        *
#*                                                                       *
#*     BigClust software is covered by the MIT license.                  *
#*                                                                       *
#*************************************************************************
# 
FF = gfortran -O3
OPEN = 
OBJ = parameters.o initisic.o clusteringmod.o objfun.o subpro.o islmb.o isic.o

all:	isic


parameters.o: parameters.f95
	$(FF) -c $(OPEN) parameters.f95

initclust.mod: r_precision.mod initisic.o initisic.f95
	$(FF) -c $(OPEN) initisic.f95

initislmb.mod: r_precision.mod initclust.mod initisic.o initisic.f95
	$(FF) -c $(OPEN) initisic.f95

initisic.o: r_precision.mod initisic.f95
	$(FF) -c $(OPEN) initisic.f95

clusteringmod.o: r_precision.mod param.mod initclust.mod clusteringmod.f95
	$(FF) -c $(OPEN) clusteringmod.f95

objfun.o: r_precision.mod initclust.mod objfun.f95
	$(FF) -c $(OPEN) objfun.f95

subpro.o: r_precision.mod param.mod subpro.f95
	$(FF) -c $(OPEN) subpro.f95

islmb.o: r_precision.mod param.mod initclust.mod initislmb.mod exe_time.mod obj_fun.mod subpro.mod islmb.f95
	$(FF) -c $(OPEN) islmb.f95

isic.o: r_precision.mod param.mod exe_time.mod initclust.mod initislmb.mod clusteringmod.mod islmb_mod.mod isic.f95
	$(FF) -c $(OPEN) isic.f95

isic: $(OBJ)
	$(FF) -o isic $(OPEN) $(OBJ)

clean:	
	rm isic $(OBJ) *.mod
	echo Clean done	