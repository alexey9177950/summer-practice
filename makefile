default:
	make get_data
	make get_iter_data
get_data:
	clang++ -std=c++14 test.cpp -o get_data
	./get_data
get_iter_data:
	clang++ -std=c++14 test_iter.cpp -o get_iter_data
	./get_iter_data
