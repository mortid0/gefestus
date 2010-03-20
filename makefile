source_dirs := . src 

search_wildcards := $(addsuffix /*.c,$(source_dirs)) 

main: $(notdir $(patsubst %.c,%.o,$(wildcard $(search_wildcards))))
	cc $^ -o $@ -lm -lc -O3
	rm *.d
	rm *.o


VPATH := $(source_dirs)
     
%.o: %.c
	cc -c -g -MD -O3 -Wall $(addprefix -I,$(source_dirs)) $< 

clean:
	rm *.d
	rm *.o

include $(wildcard *.d) 
