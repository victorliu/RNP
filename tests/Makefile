CXX = g++
CFLAGS = -O0 -ggdb

TESTPROGS = \
	test_jdqz \
	test_ira_sort \
	test_ira_zndrv1 \
	test_ira_zndrv2 \
	test_ira_zndrv3 \
	test_ira_zndrv4 \
	test_idrs \
	test_bicgstabl \
	test_gmres \
	test_lsqr \
	test_zggev \
	test_zgges \
	test_zgeev \
	test_dense_linsolve \
	test_random \
	test_io

# lsqr jdqz ira idrs gmres bicgstab
all: $(TESTPROGS)

test_io: test_io.cpp
	$(CXX) $(CFLAGS) test_io.cpp -I../inc -o test_io

test_random: test_random.cpp
	$(CXX) $(CFLAGS) test_random.cpp -I../inc -o test_random

test_dense_linsolve: test_dense_linsolve.cpp
	$(CXX) $(CFLAGS) test_dense_linsolve.cpp -I../inc -o test_dense_linsolve

test_zgeev: test_zgeev2.cpp
	$(CXX) $(CFLAGS) test_zgeev2.cpp ../src/Eigensystems.cpp -I../inc -o test_zgeev

test_zgges: test_zgges.cpp
	$(CXX) $(CFLAGS) test_zgges.cpp ../src/GeneralizedEigensystems.cpp -I../inc -o test_zgges

test_zggev: test_zggev.cpp
	$(CXX) $(CFLAGS) test_zggev.cpp ../src/GeneralizedEigensystems.cpp -I../inc -o test_zggev

test_lsqr: test_lsqr.cpp
	$(CXX) $(CFLAGS) test_lsqr.cpp -I../inc -o test_lsqr

test_gmres: test_gmres.cpp ../src/GMRES.cpp
	$(CXX) $(CFLAGS) test_gmres.cpp ../src/GMRES.cpp -I../inc -o test_gmres

test_bicgstabl: test_bicgstabl.cpp ../src/BiCGSTABl.cpp
	$(CXX) $(CFLAGS) test_bicgstabl.cpp ../src/BiCGSTABl.cpp -I../inc -o test_bicgstabl

test_idrs: test_idrs.cpp ../src/IDRs.cpp
	$(CXX) $(CFLAGS) test_idrs.cpp ../src/IDRs.cpp -I../inc -o test_idrs

test_ira_zndrv1: test_ira_zndrv1.cpp ../src/IRA.cpp
	$(CXX) $(CFLAGS) test_ira_zndrv1.cpp ../src/IRA.cpp ../src/Eigensystems.cpp -I../inc -o test_ira_zndrv1

test_ira_zndrv2: test_ira_zndrv2.cpp ../src/IRA.cpp
	$(CXX) $(CFLAGS) test_ira_zndrv2.cpp ../src/IRA.cpp ../src/Eigensystems.cpp -I../inc -o test_ira_zndrv2

test_ira_zndrv3: test_ira_zndrv3.cpp ../src/IRA.cpp
	$(CXX) $(CFLAGS) test_ira_zndrv3.cpp ../src/IRA.cpp ../src/Eigensystems.cpp -I../inc -o test_ira_zndrv3

test_ira_zndrv4: test_ira_zndrv4.cpp ../src/IRA.cpp
	$(CXX) $(CFLAGS) test_ira_zndrv4.cpp ../src/IRA.cpp ../src/Eigensystems.cpp -I../inc -o test_ira_zndrv4

test_ira_sort: test_ira_sort.cpp ../src/IRA.cpp
	$(CXX) $(CFLAGS) test_ira_sort.cpp ../src/IRA.cpp ../src/Eigensystems.cpp -I../inc -o $@

test_jdqz: test_jdqz.cpp ../src/JDQZ.cpp ../src/GeneralizedEigensystems.cpp ../src/GMRES.cpp ../src/IDRs.cpp
	$(CXX) $(CFLAGS) test_jdqz.cpp ../src/JDQZ.cpp ../src/GeneralizedEigensystems.cpp ../src/GMRES.cpp ../src/IDRs.cpp -I../inc -o test_jdqz

clean:
	rm -f $(TESTPROGS)
