/*
 * Demonstrate external client inputing and receiving outputs from a SPDZ process,
 * following the protocol described in https://eprint.iacr.org/2015/1006.pdf.
 *
 * Provides a client to bankers_bonus.mpc program to calculate which banker pays for lunch based on
 * the private value annual bonus. Up to 8 clients can connect to the SPDZ engines running
 * the bankers_bonus.mpc program.
 *
 * Each connecting client:
 * - sends a unique id to identify the client
 * - sends an integer input (bonus value to compare)
 * - sends an integer (0 meaining more players will join this round or 1 meaning stop the round and calc the result).
 *
 * The result is returned authenticated with a share of a random value:
 * - share of winning unique id [y]
 * - share of random value [r]
 * - share of winning unique id * random value [w]
 *   winning unique id is valid if ∑ [y] * ∑ [r] = ∑ [w]
 *
 * No communications security is used.
 *
 * To run with 2 parties / SPDZ engines:
 *   ./Scripts/setup-online.sh to create triple shares for each party (spdz engine).
 *   ./compile.py bankers_bonus
 *   ./Scripts/run-online bankers_bonus to run the engines.
 *
 *   ./bankers-bonus-client.x 123 2 100 0
 *   ./bankers-bonus-client.x 456 2 200 0
 *   ./bankers-bonus-client.x 789 2 50 1
 *
 *   Expect winner to be second client with id 456.
 */

#include "Math/gfp.h"
#include "Math/gf2n.h"
#include "Networking/sockets.h"
#include "Tools/int.h"
#include "Math/Setup.h"
#include "Protocols/fake-stuff.h"
#include <cstdarg>
#include <ctime>
#include <cstdlib>
#include <sodium.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort

//
// Created by wuyuncheng on 20/11/19.
//

#include "math.h"

#define SPDZ_FIXED_PRECISION 8
#define MAX_SPLIT_NUM 8
#define SPLIT_PERCENTAGE 0.8
#define LOGGER_HOME "/home/wuyuncheng/Documents/projects/CollaborativeML/log/"
FILE * logger_out;

std::string get_timestamp_str() {

    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,sizeof(buffer),"%d%m%Y%H%M%S",timeinfo);
    std::string str(buffer);

    return str;
}

void logger(FILE *out, const char *format, ...) {

    char buf[BUFSIZ] = {'\0'};
    char date_buf[50] = {'\0'};

    va_list ap;
    va_start(ap, format);
    vsprintf(buf, format, ap);
    va_end(ap);

    time_t current_time;
    current_time = time(NULL);
    struct tm *tm_struct = localtime(&current_time);
    sprintf(date_buf,"%04d-%02d-%02d %02d:%02d:%02d",
            tm_struct->tm_year + 1900,
            tm_struct->tm_mon + 1,
            tm_struct->tm_mday,
            tm_struct->tm_hour,
            tm_struct->tm_min,
            tm_struct->tm_sec);

    fprintf(out, "%s %s", date_buf, buf);
    fflush(out);
}

std::vector<int> setup_sockets(int n_parties, int my_client_id, std::vector<std::string> host_names, int port_base) {

    // Setup connections from this client to each party socket
    std::vector<int> sockets(n_parties);
    for (int i = 0; i < n_parties; i++)
    {
        set_up_client_socket(sockets[i], host_names[i].c_str(), port_base + i);
        send(sockets[i], (octet*) &my_client_id, sizeof(int));
//        octetStream os;
//        os.store(finish);
//        os.Send(sockets[i]);
        cout << "set up for " << i << "-th party succeed" << ", sockets = " << sockets[i] << ", port_num = " << port_base + i << endl;
    }
    cout << "Finish setup socket connections to SPDZ engines." << endl;
    return sockets;
}


void send_public_parameters(int type, int global_split_num, int classes_num, std::vector<int>& sockets, int n_parties) {

    octetStream os;

    gfp x = type;
    gfp y = global_split_num;
    gfp z = classes_num;

    x.pack(os);
    y.pack(os);
    z.pack(os);

    for (int i = 0; i < n_parties; i++) {
        os.Send(sockets[i]);
    }
}

// for int values
void send_private_inputs(const std::vector<gfp>& values, std::vector<int>& sockets, int n_parties)
{
    int num_inputs = values.size();
    octetStream os;
    std::vector< std::vector<gfp> > triples(num_inputs, vector<gfp>(3));
    std::vector<gfp> triple_shares(3);

    // Receive num_inputs triples from SPDZ
    for (int j = 0; j < n_parties; j++)
    {
        os.reset_write_head();
        os.Receive(sockets[j]);

        for (int j = 0; j < num_inputs; j++)
        {
            for (int k = 0; k < 3; k++)
            {
                triple_shares[k].unpack(os);
                triples[j][k] += triple_shares[k];
            }
        }
    }

    // Check triple relations (is a party cheating?)
    for (int i = 0; i < num_inputs; i++)
    {
        if (triples[i][0] * triples[i][1] != triples[i][2])
        {
            cerr << "Incorrect triple at " << i << ", aborting\n";
            exit(1);
        }
    }

    // Send inputs + triple[0], so SPDZ can compute shares of each value
    os.reset_write_head();
    for (int i = 0; i < num_inputs; i++)
    {
        gfp y = values[i] + triples[i][0];
        y.pack(os);
    }
    for (int j = 0; j < n_parties; j++)
        os.Send(sockets[j]);
}

// for float values
void send_private_batch_shares(std::vector<float> shares, std::vector<int>& sockets, int n_parties) {

    int number_inputs = shares.size();
    std::vector<long> long_shares(number_inputs);

    // step 1: convert to int or long according to the fixed precision
    for (int i = 0; i < number_inputs; ++i) {
        long_shares[i] = static_cast<int>(round(shares[i] * pow(2, SPDZ_FIXED_PRECISION)));
        //cout << "long_shares[i] = " << long_shares[i] << endl;
    }

    // step 2: convert to the gfp value and call send_private_inputs
    // Map inputs into gfp
    vector<gfp> input_values_gfp(number_inputs);
    for (int i = 0; i < number_inputs; i++) {
        input_values_gfp[i].assign(long_shares[i]);
    }

    // Run the computation
    send_private_inputs(input_values_gfp, sockets, n_parties);
    // cout << "Sent private inputs to each SPDZ engine, waiting for result..." << endl;
}

void initialise_fields(const string& dir_prefix)
{
    int lg2;
    bigint p;

    string filename = dir_prefix + "Params-Data";
    cout << "loading params from: " << filename << endl;

    ifstream inpf(filename.c_str());
    if (inpf.fail()) { throw file_error(filename.c_str()); }
    inpf >> p;
    inpf >> lg2;

    inpf.close();

    gfp::init_field(p);
    gf2n::init_field(lg2);
}

// deprecated
int receive_index(std::vector<int>& sockets)
{
    cout << "Receive best split index from the SPDZ engines" <<endl;

    octetStream os;
    os.reset_write_head();
    os.Receive(sockets[0]);

    gfp aa;
    aa.unpack(os);
    bigint index;
    to_signed_bigint(index, aa);
    int best_split_index = index.get_si();

    return best_split_index;
}


std::vector<float> receive_result(std::vector<int>& sockets, int n_parties, int size, int & best_split_index)
{
    cout << "Receive result from the SPDZ engine" << endl;
    std::vector<gfp> output_values(size);
    octetStream os;
    for (int i = 0; i < n_parties; i++)
    {
        os.reset_write_head();
        os.Receive(sockets[i]);
        for (int j = 0; j < size; j++)
        {
            gfp value;
            value.unpack(os);
            output_values[j] += value;
        }
    }

    std::vector<float> res_shares(size - 1);

    for (int i = 0; i < size - 1; i++) {
        gfp val = output_values[i];
        bigint aa;
        to_signed_bigint(aa, val);
        long t = aa.get_si();
        //cout<< "i = " << i << ", t = " << t <<endl;
        res_shares[i] = static_cast<float>(t * pow(2, -SPDZ_FIXED_PRECISION));
    }

    gfp index = output_values[size - 1];
    bigint index_aa;
    to_signed_bigint(index_aa, index);
    best_split_index = index_aa.get_si();

    return res_shares;
}


std::vector< std::vector<float> > read_training_data(int client_id, int & feature_num, int & sample_num, std::string data) {

    std::string s1("/home/wuyuncheng/Documents/projects/CollaborativeML/data/");
    std::string s11 = s1 + data + "/";
    std::string s2 = std::to_string(client_id);
    std::string data_file = s11 + "client_" + s2 + ".txt";

    std::vector <std::vector<float>> local_data;

    // read local dataset
    std::ifstream data_infile(data_file);

    if (!data_infile) {
        cout << "open %s error\n";
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(data_infile, line)) {
        std::vector<float> items;
        std::istringstream ss(line);
        std::string item;
        // split line with delimiter, default ','
        while (getline(ss, item, ',')) {
            items.push_back(::atof(item.c_str()));
        }
        local_data.push_back(items);
    }
    sample_num = local_data.size();
    feature_num = local_data[0].size();

    return local_data;
}


std::vector<int> sort_indexes(const std::vector<float> &v) {

    // initialize original index locations
    std::vector<int> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // sort 100000 running time 30ms, sort 10000 running time 3ms
    sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}


std::vector<float> compute_distinct_values(std::vector<float> feature_values, std::vector<int> sorted_indexes) {

    // now the feature values are sorted, the sorted indexes are stored in sorted_indexes
    int sample_num = feature_values.size();
    int distinct_value_num = 0;
    std::vector<float> distinct_values;
    for (int i = 0; i < sample_num; i++) {
        // if no value has been added, directly add a value
        if (distinct_value_num == 0) {
            distinct_values.push_back(feature_values[sorted_indexes[i]]);
            distinct_value_num++;
        } else {
            if (distinct_values[distinct_value_num-1] == feature_values[sorted_indexes[i]]) {
                continue;
            } else {
                distinct_values.push_back(feature_values[sorted_indexes[i]]);
                distinct_value_num++;
            }
        }
    }
    return distinct_values;
}


std::vector<float> compute_splits(std::vector<float> feature_values) {


    // use quantile sketch method, that computes k split values such that
    // there are k + 1 bins, and each bin has almost same number samples

    // basically, after sorting the feature values, we compute the size of
    // each bin, i.e., n_sample_per_bin = n/(k+1), and samples[0:n_sample_per_bin]
    // is the first bin, while (value[n_sample_per_bin] + value[n_sample_per_bin+1])/2
    // is the first split value, etc.

    // note: currently assume that the feature values is sorted, treat categorical
    // feature as label encoder sortable values

    std::vector<float> split_params;

    std::vector<int> sorted_indexes = sort_indexes(feature_values);
    int n_samples = feature_values.size();
    std::vector<float> distinct_values = compute_distinct_values(feature_values, sorted_indexes);

    // if distinct values is larger than max_bins + 1, treat as continuous feature
    // otherwise, treat as categorical feature
    if (distinct_values.size() >= MAX_SPLIT_NUM + 1) {
        // treat as continuous feature, find splits using quantile method (might not accurate when the values are imbalanced)
        split_params.push_back(MAX_SPLIT_NUM);
        int n_sample_per_bin = n_samples / (MAX_SPLIT_NUM + 1);
        for (int i = 0; i < MAX_SPLIT_NUM; i++) {
            float split_value_i = (feature_values[sorted_indexes[(i + 1) * n_sample_per_bin]]
                                   + feature_values[sorted_indexes[(i + 1) * n_sample_per_bin + 1]])/2;
            split_params.push_back(split_value_i);
        }
    }
    else if (distinct_values.size() > 1) {
        // the split values are same as the distinct values
        split_params.push_back(distinct_values.size() - 1);
        for (int i = 0; i < MAX_SPLIT_NUM; i++) {
            if (i <= (int) distinct_values.size() - 1) {
                split_params.push_back(distinct_values[i]);
            } else {
                split_params.push_back(-1);
            }
        }
    }
    else {
        // the distinct values is equal to 1, which is suspicious for the input dataset
        cout << "This feature has only one distinct value, please check it again\n";
        split_params.push_back(0);
        for (int i = 0; i < MAX_SPLIT_NUM; i++) {
//            if (i < 1) {
//                split_params.push_back(distinct_values[0]);
//            } else {
//                split_params.push_back(-1);
//            }
            split_params.push_back(-1);
        }
    }
    return split_params;
}


int main(int argc, char** argv)
{
    int my_client_id;
    int nparties;
    int port_base = 20000;
    std::vector<std::string> host_names;
    //string host_name = "localhost";
    host_names.push_back("127.0.0.1");
    host_names.push_back("127.0.0.1");
    host_names.push_back("127.0.0.1");

    if (argc < 3) {
        cout << "Usage is bankers-bonus-client <client identifier> <number of spdz parties> "
           << "<salary to compare> <finish (0 false, 1 true)> <optional host name, default localhost> "
           << "<optional spdz party port base number, default 14000>" << endl;
        exit(0);
    }

    my_client_id = atoi(argv[1]);
    nparties = atoi(argv[2]);
    std::string data_file = "m3_n10000_d30_c6";
    if (argc > 3)
        data_file = argv[3];
    if (argc > 4)
        port_base = atoi(argv[4]);

    std::string logger_file_name = LOGGER_HOME;

    logger_file_name += data_file;
    logger_file_name += "_";
    logger_file_name += get_timestamp_str();
    logger_file_name += "_client";
    logger_file_name += to_string(my_client_id);
    logger_file_name += ".txt";
    logger_out = fopen(logger_file_name.c_str(), "wb");


    struct timeval spdz_training_1, spdz_training_2;
    double spdz_training_time = 0;
    gettimeofday(&spdz_training_1, NULL);

    // init static gfp
    string prep_data_prefix = get_prep_dir(nparties, 128, gf2n::default_degree());
    initialise_fields(prep_data_prefix);
    bigint::init_thread();

    cout << "Begin setup sockets" << endl;

    // Setup connections from this client to each party socket
    vector<int> sockets = setup_sockets(nparties, my_client_id, host_names, port_base);
    cout << "sockets[0] = " << sockets[0] << endl;
    cout << "sockets[1] = " << sockets[1] << endl;
    cout << "sockets[2] = " << sockets[2] << endl;

    cout << "Finish setup socket connections to SPDZ engines." << endl;

    int feature_num, sample_num;
    std::vector< std::vector<float> > local_data = read_training_data(my_client_id, feature_num, sample_num, data_file);

    cout << "Correct read training data" << endl;
    cout << "sample_num = " << sample_num << endl;
    cout << "feature_num = " << feature_num << endl;

    // only use the former 80%
    std::vector< std::vector<float> > training_data;

    if (my_client_id == 0) {

        std::vector<float> training_labels;
        int training_data_num = sample_num * SPLIT_PERCENTAGE;
        cout << "training_data_num = " << training_data_num << endl;
        for (int i = 0; i < training_data_num; i++) {
            training_labels.push_back(local_data[i][feature_num-1]);
            std::vector<float> tmp;
            for (int j = 0; j < feature_num - 1; j++) {
                tmp.push_back(local_data[i][j]);
            }
            training_data.push_back(tmp);
        }

        feature_num = feature_num - 1;

        for (int i = 0; i < (int) training_data.size(); i++) {
            for (int j = 0; j < (int) training_data[0].size(); j++) {
                std::vector<float> x;
                x.push_back(training_data[i][j]);
                send_private_batch_shares(x, sockets, nparties);
            }
        }

        cout << "Finish send training data to SPDZ engines." << endl;

        for (int i = 0; i < (int) training_labels.size(); i++) {
            std::vector<float> x;
            x.push_back(training_labels[i]);
            send_private_batch_shares(x, sockets, nparties);
        }

        cout << "Finish send training labels to SPDZ engines." << endl;

    } else {
        int training_data_num = sample_num * SPLIT_PERCENTAGE;
        cout << "training_data_num = " << training_data_num << endl;
        for (int i = 0; i < training_data_num; i++) {
            std::vector<float> tmp;
            for (int j = 0; j < feature_num; j++) {
                tmp.push_back(local_data[i][j]);
                std::vector<float> x;
                x.push_back(local_data[i][j]);
                send_private_batch_shares(x, sockets, nparties);
            }
            training_data.push_back(tmp);
        }

        cout << "Finish send training data to SPDZ engines." << endl;
    }

    // send split parameters
    for (int j = 0; j < feature_num; j++) {

        std::vector<float> feature_values;
        for (int i = 0; i < (int) training_data.size(); i++) {
            feature_values.push_back(training_data[i][j]);
        }

        std::vector<float> split_params = compute_splits(feature_values);

        if (split_params.size() != MAX_SPLIT_NUM + 1) {
            cout << "Error split params size." << endl;
        }

        for (int s = 0; s < MAX_SPLIT_NUM + 1; s++) {
            // the first one is the real split number
            // the rest are the split values
            std::vector<float> x;
            //cout << "feature " << j << ", split param　" << s << " = " << split_params[s] << endl;
            x.push_back(split_params[s]);
            send_private_batch_shares(x, sockets, nparties);
        }
    }

    cout << "Finish send split parameters to SPDZ engines." << endl;

    int finished = receive_index(sockets);

    cout << "finished = " << finished << endl;

    for (unsigned int i = 0; i < sockets.size(); i++) {
        close_client_socket(sockets[i]);
    }

    gettimeofday(&spdz_training_2, NULL);
    spdz_training_time += (double)((spdz_training_2.tv_sec - spdz_training_1.tv_sec) * 1000 +
                                   (double)(spdz_training_2.tv_usec - spdz_training_1.tv_usec) / 1000);
    logger(logger_out, "*********************************************************************");
    logger(logger_out, "******** SPDZ training time: %'.3f ms **********\n", spdz_training_time);
    logger(logger_out, "*********************************************************************");

    return 0;
}
