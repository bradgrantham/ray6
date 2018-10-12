ray6: brads.o distrib.o grf_color.o grf_math.o grf_obj.o ray6.o

lines4010: lines4010.o grf_math.o grf_obj.o
	$(CC) $< -o $@ -l4010

brads.o: brads.h brads.c

distrib.o: distrib.c grf_math.h

grf_color.o: grf_color.c grf_color.h

grf_math.o: grf_math.c grf_math.h brads.h

grf_obj.o: grf_obj.c grf_obj.h

ray6.o: ray6.c grf_math.h grf_view.h grf_obj.h


