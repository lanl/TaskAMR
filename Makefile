test:
	./test_linear.py
	./test_euler.py

prof:
	../../github/legion/tools/legion_prof.py -o ./prof prof0

spy:  
	../../github/legion/tools/legion_spy.py -dez spy0

clean:  
	rm -r prof* spy* *txt

