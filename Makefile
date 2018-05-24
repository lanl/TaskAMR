LEGION_ROOT ?= ../../github/legion

test:
	$(LEGION_ROOT)/language/regent.py unit_tests.rg
	./test_linear_amr.py
	./test_linear.py
	./test_euler.py

prof:
	$(LEGION_ROOT)/tools/legion_prof.py -o ./prof prof0

spy:  
	$(LEGION_ROOT)/tools/legion_spy.py -dez spy0

clean:  
	rm -r prof* spy* *txt *~

