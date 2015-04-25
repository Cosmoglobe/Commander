#include <cstdlib>
#include <cstdio>
#include <unistd.h>
#include <stdint.h>
#include <signal.h>
#include <cstring>
#include <string>
#include <cstdarg>
#include <map>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <pthread.h>
#include <fenv.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
using namespace std;

typedef int (*forthandler)(int*);

map<int,forthandler> handlers;
void sig_wrapper(int signal)
{
	forthandler h = handlers[signal];
	if(!h)
	{
		fprintf(stderr, "Caught signal %d, but don't know how to handle it!\n", signal);
		exit(1);
	}
	int survive = h(&signal);
	if(!survive) exit(1);
}
int nfiles = 0;
map<int,FILE*> files;

extern "C"
{
	int get_pid_() { return getpid(); }
	int64_t get_mem_use_()
	{
		stringstream path; path << "/proc/" << getpid() << "/statm";
		FILE * f = fopen(path.str().c_str(),"r");
		int64_t npage;
		fscanf(f, "%lld", &npage);
		fclose(f);
		return npage*getpagesize();
	}
	int64_t get_max_mem_use_()
	{
		stringstream path; path << "/proc/" << getpid() << "/status";
		ifstream f(path.str().c_str());
		string line;
		int64_t vm_peak_kb;
		while(getline(f, line))
			if(sscanf(line.c_str(), "VmPeak: %ld kB\n", &vm_peak_kb) == 1)
				return vm_peak_kb << 10;
		return 0;
	}
	void set_sig_handler_(int * signal, forthandler foo) {
		handlers[*signal] = foo;
		struct sigaction ho, h;
		h.sa_handler = sig_wrapper;
		h.sa_flags = 0;
		sigemptyset(&h.sa_mask);
		sigaction(*signal, &h, &ho);
	}
	void clear_sig_handler_(int * signal) {
		struct sigaction ho, h;
		handlers[*signal] = 0;
		h.sa_handler = SIG_DFL;
		h.sa_flags = 0;
		sigemptyset(&h.sa_mask);
		sigaction(*signal, &h, &ho);
	}
	void wall_time_(double * t) {
		timeval tv;
		gettimeofday(&tv,0);
		*t = tv.tv_sec + 1e-6*tv.tv_usec;
	}
	void usleep_(int * t) { usleep(*t); }
	// Random-access-oriented file that supports locking over NFS.
	// Probably slow, so do not use for large amounts of data.
	int open_atomic_file_(char * fname, int * trunc, int len)
	{
		string filename(fname,fname+len);
		int mode = O_RDWR | O_CREAT | O_SYNC;
		if(*trunc) mode |= O_TRUNC;
		return open(filename.c_str(), mode, 0666);
	}
	void close_atomic_file_(int * fd) { close(*fd); }
	void lock_atomic_file_(int * fd, int * pos, int * len)
	{
		lseek(*fd, *pos, SEEK_SET);
		lockf(*fd, F_LOCK, *len);
	}
	void unlock_atomic_file_(int * fd, int * pos, int * len)
	{
		lseek(*fd, *pos, SEEK_SET);
		lockf(*fd, F_ULOCK, *len);
	}
	void read_atomic_file_(int * fd, int * pos, int * len, void * data)
	{
		lseek(*fd, *pos, SEEK_SET);
		read(*fd, data, *len);
	}
	void write_atomic_file_(int * fd, int * pos, int * len, void * data)
	{
		lseek(*fd, *pos, SEEK_SET);
		write(*fd, data, *len);
	}
	void append_atomic_file_(int * fd, int * len, void * data)
	{
		lseek(*fd, 0, SEEK_END);
		write(*fd, data, *len);
	}
	void thread_create_(void *(*fun)(void*), void * args, pthread_t * thread, int * status)
	{
		*status = pthread_create(thread, NULL, fun, &args);
	}
	void thread_exit_() { pthread_exit(NULL); }
	void thread_join_(pthread_t * thread) { pthread_join(*thread, NULL); }
	int nfork_(int * n)
	{
		int i; for(i = 1; i < *n && fork(); i++);
		return i-1;
	}
	int popen_(const char * command, const char * type, int len1, int len2)
	{
		FILE * f = popen(string(command,command+len1).c_str(),string(type,type+len2).c_str());
		if(!f) return 0;
		nfiles++;
		files[nfiles] = f;
		return nfiles;
	}
	void pclose_(int * i) {
		map<int,FILE*>::iterator j = files.find(*i);
		if(j == files.end()) return;
		fclose(j->second);
		files.erase(j);
	}
	void preadline_(int * i, char * buffer, int * n, int buflen)
	{
		FILE * f = files[*i];
		if(!fgets(buffer, buflen, f)) *n = -1;
		else
		{
			*n = strlen(buffer);
			if(*n > 1 && buffer[*n-1] == '\n') buffer[*n-1] = ' ';
			for(int j = *n; j < buflen; j++)
				buffer[j] = ' ';
		}
	}
	void pwriteline_(int * i, char * buffer, int buflen)
	{
		FILE * f = files[*i];
		string s(buffer,buffer+buflen);
		s += '\n';
		fputs(s.c_str(), f);
	}
	void getenv_(const char * name, char * value, int l1, int l2)
	{
		char * s = getenv(string(name,name+l1).c_str());
		int n = 0;
		if(s) {
			n = strlen(s);
			strncpy(value, s, l2);
		}
		for(int i = n; i < l2; i++) value[i] = ' ';
	}
	void setenv_(const char * name, char * value, int l1, int l2)
	{
		setenv(string(name,name+l1).c_str(), string(value,value+l2).c_str(),1);
	}
	void get_hostname_(char * value, int l2)
	{
		utsname u; uname(&u);
		int n = 0;
		n = strlen(u.nodename);
		strncpy(value, u.nodename, l2);
		for(int i = n; i < l2; i++) value[i] = ' ';
	}
	void truncate_(const char * name, int * len, int l)
	{
		truncate(string(name,name+l).c_str(), *len);
	}
	void toupper_(char * str, int l)
	{
		for(int i = 0; i < l; i++) str[i] = toupper(str[i]);
	}
	void tolower_(char * str, int l)
	{
		for(int i = 0; i < l; i++) str[i] = tolower(str[i]);
	}
	bool is_dir_(char * path, int len)
	{
		string s(path, path+len);
		struct stat buf;
		if(stat(s.c_str(), &buf) != 0) return false;
		return S_ISDIR(buf.st_mode);
	}
	bool mkdir_c_(char * path, int len)
	{
		string s(path, path+len);
		return mkdir(s.c_str(), 0777) == 0;
	}
	bool mv_c_(char * a, char * b, int la, int lb)
	{
		return rename(string(a,a+la).c_str(), string(b,b+lb).c_str()) == 0;
	}
	bool rm_c_(char * path, int len)
	{
		return unlink(string(path,path+len).c_str()) == 0;
	}
	int ishift_(int * val, int * steps) { return *val << *steps; }

  //int svn_revision_() { return REVISION; }
        double nan = numeric_limits<double>::quiet_NaN(),
		snan = numeric_limits<double>::signaling_NaN(),
		infinity = numeric_limits<double>::infinity();
	int64_t fe_divbyzero = FE_DIVBYZERO, fe_inexact = FE_INEXACT, fe_nan = FE_INVALID, fe_overflow = FE_OVERFLOW, fe_underflow = FE_UNDERFLOW;
	void fe_enable_(int64_t * a) { feenableexcept(*a); }
	void fe_disable_(int64_t * a) { fedisableexcept(*a); }
	uint32_t count_set_bits_(uint32_t * j)
	{
		uint32_t i = *j;
		i = i - ((i >> 1) & 0x55555555);
		i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
		return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
	}
}
