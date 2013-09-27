#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cstring>

#include <sstream>
#include <vector>
#include <string>

#include <stdint.h>
#include <sys/stat.h>
#include <functional>

#include <getopt.h>

#include "Common.hpp"
#include "CountBF.hpp"

#include "HashTables.hpp"
#include "fastq.hpp"
#include "Kmer.hpp"
#include "KmerIterator.hpp"
#include "BloomFilter.hpp"
#include "bloom_filter.hpp"

#include "QLogTable.hpp"

// structs for getopt
struct CountBF_ProgramOptions {
  size_t k;
  size_t nkmers;
  string output;
  bool verbose;
  bool quake;
  size_t bf;
  size_t qs;
  size_t threads;
  uint32_t seed;
  size_t read_chunksize;
  vector<string> files;

  CountBF_ProgramOptions() : k(0), nkmers(0), verbose(false), quake(false), bf(4), qs(0), threads(1), seed(0), read_chunksize (10000) {}
};

void CountBF_PrintUsage() {
  cerr << "BFCounter " << BFC_VERSION << endl << endl;
  cerr << "Counts occurrences of k-mers in fastq or fasta files and saves results" << endl << endl;
  cerr << "Usage: BFCounter count [options] ... FASTQ files";
  cerr << endl << endl <<
    "-k, --kmer-size=INT             Size of k-mers, at most " << (int) (Kmer::MAX_K-1)<< endl << 
    "-n, --num-kmers=LONG            Estimated number of k-mers (upper bound)" << endl <<
    "-t, --threads=INT               Number of threads to use (default 1)" << endl <<
    "-c, --chunk-size=INT            Number of reads to proccess in parallel (default 10000)" << endl << 
    "-s, --seed=INT                  seed (32-bit int) for randomization" << endl <<
    "-o, --output=STRING             Filename for output" << endl <<
    "-b, --bloom-bits=INT            Number of bits to use in Bloom filter (default=4)" << endl <<
    "    --quake                     Count q-mers for use with Quake (default=FALSE)" << endl <<
    "    --quality-scale=INT         Quality-scale used in Quake mode, only 64 and 33 (default=64)" << endl <<
    "    --verbose                   Print lots of messages during run" << endl << endl
    ;
}




void CountBF_ParseOptions(int argc, char **argv, CountBF_ProgramOptions &opt) {
  int verbose_flag = 0;
  int quake_flag = 0;
  const char* opt_string = "n:k:o:b:t:s:c:";
  static struct option long_options[] =
  {
    {"verbose", no_argument,  &verbose_flag, 1},
    {"kmer-size", required_argument, 0, 'k'},
    {"num-kmers", required_argument, 0, 'n'},
    {"seed", required_argument, 0, 's'},
    {"threads", required_argument, 0, 't'},
    {"chunk-size", required_argument, 0, 'c'},
    {"output", required_argument, 0, 'o'},
    {"bloom-bits", required_argument, 0, 'b'},
    {"quality-scale", required_argument, 0, 0},
    {"quake", no_argument, &quake_flag, 1},
    {0,0,0,0}
  };

  int option_index = 0;
  int c;
  stringstream ss;
  while (true) {
    c = getopt_long(argc,argv,opt_string, long_options, &option_index);

    if (c == -1) {
      break;
    }

    switch (c) {
    case 0: 
      if (strcmp("quality-scale", long_options[option_index].name) == 0) {
	opt.qs = atoi(optarg);
      }
      break;
    case 'k': 
      opt.k = atoi(optarg); 
      break;
    case 'o': 
      opt.output = optarg;
      break;
    case 'n': 
      ss << optarg;
      ss >> opt.nkmers;
      break;
    case 'b':
      opt.bf = atoi(optarg);
      break;
    case 's':
      opt.seed = atoi(optarg);
      break;
    case 't':
      opt.threads = atoi(optarg);
      break;
    case 'c':
      opt.read_chunksize = atoi(optarg);
      break;
    default: break;
    }
  }

  // all other arguments are fast[a/q] files to be read
  for (int i = optind; i < argc; i++) {
    opt.files.push_back(argv[i]);
  }
  
  if (verbose_flag) {
    opt.verbose = true;
  }
  if (quake_flag) {
    opt.quake = true;
  }
}


bool GuessQualityScore(CountBF_ProgramOptions &opt) {
  opt.qs = 64;
  
  FastqFile FQ(opt.files);
  char name[8196], s[8196], qual[8196];
  size_t name_len, len;
  size_t nread = 0;
  // for each read

  uint8_t min = 255, max = 0;
  while (FQ.read_next(name, &name_len, s, &len, NULL, qual) >= 0) {
    nread++;
    // check all quality scores 

    for (size_t i = 0; i < len; ++i) {
      min = (qual[i] < min) ? qual[i] : min;
      max = (qual[i] > max) ? qual[i] : max;
    }
    if (nread > 10000) {
      break;
    }
  }
  FQ.close();

  if (min < '!' || max > '~' || max-min > 62) {
    return false;
  } 

  if (min < 64) {
    opt.qs = 33;
  }
  return true;

}

bool CountBF_CheckOptions(CountBF_ProgramOptions &opt) {
  bool ret = true;

  if (opt.k <= 0 || opt.k >= MAX_KMER_SIZE) {
    cerr << "Error, invalid value for kmer-size: " << opt.k << endl;
    cerr << "Values must be between 1 and " << (MAX_KMER_SIZE-1) << endl;
    ret = false;
  }

  if (opt.nkmers <= 0) {
    cerr << "Error, invalid value for num-kmers: " << opt.nkmers << endl;
    cerr << "Values must be positive integers" << endl;
    ret = false;
  }

  if (opt.threads <= 0) {
    cerr << "Error, invalid value for threads: " << opt.threads << endl;
    cerr << "Values must be positive integers" << endl;
    ret = false;
  }

 if (opt.read_chunksize <= 0) {
    cerr << "Error, invalid value for chunk-size: " << opt.read_chunksize << endl;
    cerr << "Values must be positive integers" << endl;
    ret = false;
  }



  if (opt.files.size() == 0) {
    cerr << "Need to specify files for input" << endl;
    ret = false;
  } else {
    struct stat stFileInfo;
    vector<string>::const_iterator it;
    int intStat;
    for(it = opt.files.begin(); it != opt.files.end(); ++it) {
      intStat = stat(it->c_str(), &stFileInfo);
      if (intStat != 0) {
	cerr << "Error: file not found, " << *it << endl;
	ret = false;
      }
    }
  }
  
  if (opt.bf <= 0) {
    cerr << "Invalid value for bloom filter size" << endl;
    ret = false;
  }

  if (opt.quake) {
    if (opt.qs != 0) {
      if (opt.qs != 64 && opt.qs != 33) {
	cerr << "Invalid value for quality-scale, we only accept 64 and 33" << endl;
	ret = false;
      }
    } else {
      // we'll guess the quality score
      if (opt.verbose) {
	cerr << "Guessing quality scale:";
      }
      if (!GuessQualityScore(opt)) {
	cerr << endl << "Could not guess quality scale from sequence reads, set manually to 33 or 64" << endl;
	ret = false;
      }	else if (opt.verbose) {
	cerr << " quality scale " << opt.qs << endl;
      }
    }
  }

  //TODO: check if we have permission to write to outputfile
  
  return ret;

}

void CountBF_PrintSummary(const CountBF_ProgramOptions &opt) {
  cerr << "Using bloom filter size: " << opt.bf << " bits" << endl;
  cerr << "Estimated false positive rate: ";
  double fp = pow(pow(.5,log(2.0)),(double) opt.bf);
  cerr << fp << endl;
}

void CountBF_Quake(const CountBF_ProgramOptions &opt) {
  // create hash table and bloom filter
  
  size_t k = Kmer::k;
  hmapq_t kmap;
  bloom_filter BF(opt.nkmers, (size_t) opt.bf, (unsigned long) time(NULL));

  char name[8196], s[8196], qual[8196];
  size_t name_len, len;
  
  uint64_t n_read = 0;
  uint64_t num_kmers = 0;
  uint64_t filtered_kmers =0 ;
  uint64_t total_cov = 0;

  // loops over all files
  FastqFile FQ(opt.files);

  // for each read
  while (FQ.read_next(name, &name_len, s, &len, NULL, qual) >= 0) {
    if (len < k) {
      continue;
    }
    Kmer km(s);
    for (size_t i = 0; i <= len-k; ++i) {
      num_kmers++;
      if (i > 0) {
	km = km.forwardBase(s[i+k-1]);
      }
    
      Kmer tw = km.twin();
      Kmer rep = (km < tw) ? km : tw;
      if (BF.contains(rep)) {
	// has no effect if already in map
	pair<hmapq_t::iterator, bool> ref = kmap.insert(make_pair(rep, 0.0f));
      } else {
	BF.insert(rep);
      }
    }

    ++n_read;
    if (opt.verbose && n_read % 1000000 == 0) {
      cerr << "processed " << n_read << " reads" << endl;
    }
  }
  
  if (opt.verbose) {
    cerr << "re-open all files" << endl;
  }

  FQ.reopen();
  hmapq_t::iterator it;
  float qlogsum;
  while (FQ.read_next(name, &name_len, s, &len, NULL, qual) >= 0) {
    if (len < k) {
      continue;
    }
    Kmer km(s);
    qlogsum = 0.0f;
    for (size_t i = 0; i < k; ++i) {
      qlogsum += qlogtable[(uint8_t)qual[i]-opt.qs];
    }

    for (size_t i = 0; i <= len-k; ++i) {
      if (i > 0) {
	km = km.forwardBase(s[i+k-1]);
	qlogsum += qlogtable[(uint8_t)qual[i+k-1]-opt.qs] - qlogtable[(uint8_t)qual[i-1]-opt.qs];
      }
      
      Kmer tw = km.twin();
      Kmer rep = (km < tw) ? km : tw;
      
      it = kmap.find(rep);
      if (it != kmap.end()) {
	it->second += exp(qlogsum);
	total_cov += 1;
      }
    } 
  }

  FQ.close();

  // the hash map needs an invalid key to mark as deleted
  Kmer km_del;
  km_del.set_deleted();
  kmap.set_deleted_key(km_del);

  if (opt.verbose) {
    cerr << "processed " << num_kmers << " kmers in " << n_read  << " reads"<< endl;
    cerr << "found " << kmap.size() << " non-filtered kmers, kept all" << endl;
    filtered_kmers = num_kmers - total_cov;
    
    cerr << "total coverage " << total_cov << ", estimated number of kmers " << filtered_kmers << endl;
    cerr << "average coverage " << (total_cov / ((double) kmap.size())) << endl;
  }

  if (opt.verbose) {
    cerr << "Writing hash table to file " << opt.output << " .. "; cerr.flush();
    cerr << "hashtable size is " << kmap.size()/(1<<20) << "MB"  << endl;
  }
  FILE* f = fopen(opt.output.c_str(), "wb");
  if (f == NULL) {
    cerr << "Error could not write to file!" << endl;
  } else {
    // first metadata for hash table
    kmap.write_metadata(f);
    // then the actual hashtable
    kmap.write_nopointer_data(f);
    fclose(f);
    f = NULL;
  }

  if (opt.verbose) {
    cerr << " done" << endl << endl;
    cerr << " convert the file to tabular format using the command " << endl <<
        "    BFCounter dump -k " << k << " -i " << opt.output << " -o output_file " << endl;
       
  }
    
}

void CountBF_Normal(const CountBF_ProgramOptions &opt) {
  // create hash table and bloom filter

  hmap_t kmap;
  hmapL_t kmap_Large;
  
  size_t num_threads = opt.threads;
#ifdef _OPENMP
  omp_set_num_threads(num_threads);
#endif
  
  uint32_t seed = opt.seed;
  if (seed == 0) {
    seed = (uint32_t) time(NULL);
  }
  BloomFilter BF(opt.nkmers, (size_t) opt.bf, seed);
  
  bool done = false;
  
  char name[8196],s[8196];//, qual[8196];
  size_t name_len,len;

  uint64_t n_read = 0;
  uint64_t num_kmers = 0;  
  uint64_t filtered_kmers = 0;
  uint64_t total_cov = 0;
  size_t read_chunksize = opt.read_chunksize;

  // loops over all files
  FastqFile FQ(opt.files);
  string *readv = new string[read_chunksize];
  vector<Kmer> *parray = new vector<Kmer>[num_threads];
  vector<Kmer> *smallv;
  size_t round = 0;

   
  // for each batch
  while (!done) {
    size_t reads_now = 0;
    while (reads_now < read_chunksize) {
      if (FQ.read_next(name, &name_len, s, &len, NULL, NULL) >= 0) {
	readv[reads_now].assign(s);
	++n_read;
	++reads_now;
      } else {
	done = true;
	break;
      }
    }
    ++round;

#pragma omp parallel default(shared) private(smallv) shared(parray, readv, BF, reads_now) reduction(+: num_kmers, n_read)
    {
      KmerIterator iter, iterend;
      size_t threadnum = 0;
#ifdef _OPENMP
      threadnum = omp_get_thread_num();
#endif
      smallv = &parray[threadnum];

      #pragma omp for nowait
      for (size_t index = 0; index < reads_now; ++index) {
	// for each read in our batch
	const char *cstr = readv[index].c_str();
	iter = KmerIterator(cstr);
	n_read++;
	for(; iter != iterend; ++iter) {
	  // for each valid k-mer in read
	  ++num_kmers;
	  Kmer rep = iter->first.rep();
	  size_t r = BF.search(rep);
	  if (r == 0) {
	    // in bf
	    smallv->push_back(rep);
	  } else {
	    if (BF.insert(rep) == r) {
	      // inserted by us
	    } else {
	      // might have been inserted by other thread simultaneously
	      smallv->push_back(rep);
	    }
	  }
	} // done with k-mers
      } // done with read
    } // done with this batch

    // this part is serial
    for (size_t i = 0; i < num_threads; i++) {
      for (vector<Kmer>::const_iterator it = parray[i].begin(); it != parray[i].end(); ++it) {
	kmap.insert(KmerIntPair(*it,0)); // no extra effect if duplicated
      }
      parray[i].clear();
    }

    if (opt.verbose && read_chunksize > 1) {
      cerr << "processed " << n_read << " reads" << endl;
    }
  }
  
  if (opt.verbose) {
    cerr << "re-open all files" << endl;
  }
  // close all files, reopen and get accurate counts;
  FQ.reopen();


  n_read = 0; // reset counter

  done = false;
  while (!done) {
    size_t reads_now = 0;
    while (reads_now < read_chunksize) {
      if (FQ.read_next(name, &name_len, s, &len, NULL, NULL) >= 0) {
	readv[reads_now].assign(s);
	++n_read;
	++reads_now;
      } else {
	done = true;
	break;
      }
    }
    ++round;

#pragma omp parallel default(shared) private(smallv) shared(parray, readv, BF, reads_now) reduction(+: total_cov, n_read)
    {
      hmap_t::iterator it;
      KmerIterator iter, iterend;
      size_t threadnum = 0;
#ifdef _OPENMP
      threadnum = omp_get_thread_num();
#endif
      smallv = &parray[threadnum];

#pragma omp for nowait
      for (size_t index = 0; index < reads_now; ++index) {
	// for each read in our batch
	const char *cstr = readv[index].c_str();
	iter = KmerIterator(cstr);
	n_read++;
	for(; iter != iterend; ++iter) {
	  // for each valid k-mer in read
	  Kmer rep = iter->first.rep();
	  it = kmap.find(rep);
	  if (it != kmap.end()) {
            bool b = true;
            unsigned int val = it->GetVal();
            if (val < KmerIntPair::MaxVal) {
              b = it->ParallelIncrement();
            }
	    if (!b || val == KmerIntPair::MaxVal) { // ok we did not increment is so it was 255 already
	      smallv->push_back(rep); // large values, handle serially
	    }
	    total_cov += 1;
	  }
	} // done with k-mers
      } // done with read
    } // done with this batch

    // this part is serial
    for (size_t i = 0; i < num_threads; i++) {
      for (vector<Kmer>::const_iterator it = parray[i].begin(); it != parray[i].end(); ++it) {
	Kmer rep = *it;
	hmapL_t::iterator l_it = kmap_Large.find(rep);
	if (l_it == kmap_Large.end()) {
	  kmap_Large.insert(make_pair(rep,KmerIntPair::MaxVal+1));
	} else {
          l_it->second += 1;
	}
      }
      parray[i].clear();
    }

    if (opt.verbose && read_chunksize > 1) {
      cerr << "processed " << n_read << " reads" << endl;
    }
  }
  
  FQ.close();

  if (opt.verbose) {
    cerr << "closed all files" << endl;
  }

  // the hash map needs an invalid key to mark as deleted
  Kmer km_del;
  km_del.set_deleted();
  kmap.set_deleted_key(km_del);
  size_t n_del =0 ;

  for(hmap_t::iterator it = kmap.begin(); it != kmap.end(); ) {
    if (it->GetVal() <= 1) {
      hmap_t::iterator del(it);
      ++it;
      // remove k-mer that got through the bloom filter
      kmap.erase(del);
      ++n_del;
    } else {
      ++it;
    }
  }

  total_cov -= n_del;

  if (opt.verbose) {
    cerr << "processed " << num_kmers << " kmers in " << n_read  << " reads"<< endl;
    cerr << "found " << kmap.size() << " non-filtered kmers, removed " << n_del << endl;
    filtered_kmers = num_kmers - total_cov;
    
    cerr << "total coverage " << total_cov << ", estimated number of kmers " << filtered_kmers << endl;
    cerr << "average coverage " << (total_cov / ((double) kmap.size())) << endl;

  }

  if (opt.verbose) {
    cerr << "Writing hash table to file " << opt.output << " .. "; cerr.flush();
    cerr << "hashtable size is " << kmap.size()  << " k-mers" << endl;
  }
  FILE* f = fopen(opt.output.c_str(), "wb");
  if (f == NULL) {
    cerr << "Error could not write to file!" << endl;
  } else {
    // first metadata for hash table
    kmap.write_metadata(f);
    // then the actual hashtable
    kmap.write_nopointer_data(f);
    kmap_Large.write_metadata(f);
    kmap_Large.write_nopointer_data(f);
    fclose(f);
    f = NULL;
  }
  if (opt.verbose) {
    cerr << " done" << endl << endl;
    cerr << " convert the file to tabular format using the command " << endl <<
        "    BFCounter dump -k " << Kmer::k << " -i " << opt.output << " -o output_file " << endl;
  }
}

void CountBF(int argc, char **argv) {
  
  CountBF_ProgramOptions opt;
  CountBF_ParseOptions(argc,argv,opt);

  if (argc < 2) {
    CountBF_PrintUsage();
    exit(1);
  }
  
  if (!CountBF_CheckOptions(opt)) {
    CountBF_PrintUsage();
    exit(1);
  }
  
  // set static global k-value
  Kmer::set_k(opt.k);

  if (opt.verbose) {
    CountBF_PrintSummary(opt);
  }

  if (opt.quake) {
    CountBF_Quake(opt);
  } else {
    CountBF_Normal(opt);
  }
}
