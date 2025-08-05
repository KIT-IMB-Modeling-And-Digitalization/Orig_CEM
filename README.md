Install using the following command, it will find your OS automatically:

python setup.py build_ext

Then based on your OS pick one of the files in simulations directory - either 01_linux.py or 01_win.py

At the top assign the id to your simulation and then run using the following command:

python FILENAME.py

It will automatically create a directory with correlation files and run the simulation till the end and finally moves everything to a directory named results_id within simulations directory.

Attention: To run the program in linux terminal, please remain inside linux itself and do not use the mounted directory of windows like /mnt/c/... because it will crash

