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

#include <sodium.h>
#include <iostream>
#include <sstream>
#include <fstream>

//
// Created by wuyuncheng on 20/11/19.
//

#include "math.h"

#define SPDZ_FIXED_PRECISION 16


std::vector<int> setup_sockets(int n_parties, int my_client_id, const std::string host_name, int port_base) {

    // Setup connections from this client to each party socket
    std::vector<int> sockets(n_parties);
    for (int i = 0; i < n_parties; i++)
    {
        set_up_client_socket(sockets[i], host_name.c_str(), port_base + i);
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

void send_private_batch_shares(std::vector<float> shares, std::vector<int>& sockets, int n_parties) {

    int number_inputs = shares.size();
    std::vector<int64_t> long_shares(number_inputs);

    // step 1: convert to int or long according to the fixed precision
    for (int i = 0; i < number_inputs; ++i) {
        long_shares[i] = static_cast<int64_t>(round(shares[i] * pow(2, SPDZ_FIXED_PRECISION)));
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
        int64_t t = aa.get_si();
        //cout<< "i = " << i << ", t = " << t <<endl;
        res_shares[i] = static_cast<float>(t * pow(2, -SPDZ_FIXED_PRECISION));
    }

    gfp index = output_values[size - 1];
    bigint index_aa;
    to_signed_bigint(index_aa, index);
    best_split_index = index_aa.get_si();

    return res_shares;
}


int main(int argc, char** argv)
{
    int my_client_id;
    int nparties;
    int port_base = 18000;
    string host_name = "localhost";

    if (argc < 3) {
        cout << "Usage is bankers-bonus-client <client identifier> <number of spdz parties> "
           << "<salary to compare> <finish (0 false, 1 true)> <optional host name, default localhost> "
           << "<optional spdz party port base number, default 14000>" << endl;
        exit(0);
    }

    my_client_id = atoi(argv[1]);
    nparties = atoi(argv[2]);
    if (argc > 3)
        host_name = argv[3];
    if (argc > 4)
        port_base = atoi(argv[4]);

    // init static gfp
    string prep_data_prefix = get_prep_dir(nparties, 128, gf2n::default_degree());
    initialise_fields(prep_data_prefix);
    bigint::init_thread();

    cout<<"Begin setup sockets"<<endl;

    for (int xx = 0; xx < 15; xx++) {

        cout << "iteration xx = " << xx << endl;

        // Setup connections from this client to each party socket
        vector<int> sockets = setup_sockets(nparties, my_client_id, host_name, port_base);
        cout << "sockets[0] = " << sockets[0] << endl;
        cout << "sockets[1] = " << sockets[1] << endl;
        cout << "sockets[2] = " << sockets[2] << endl;


        cout << "Finish setup socket connections to SPDZ engines." << endl;

        int type = 0;
        int global_split_num = 250;
        int classes_num = 2;

        if (my_client_id == 0) {
            send_public_parameters(type, global_split_num, classes_num, sockets, nparties);
            cout << "Finish send public parameters to SPDZ engines." << endl;
        }

        cout << "Ready to send private statistics" << endl;
        vector< vector<float> > statistics;
        for (int i = 0; i < global_split_num; i++) {
            vector<float> tmp;
            for (int j = 0; j < classes_num * 2; j++) {
                tmp.push_back(my_client_id * 0.1 + (i + 1) * 0.01 + (j + 1) * 0.005);
            }
            statistics.push_back(tmp);
        }

        for (int i = 0; i < global_split_num; i++) {
            for (int j = 0; j < classes_num * 2; j++) {
                vector<float> x;
                x.push_back(statistics[i][j]);
                //cout << "statistics["<< i << "][" << j << "] = "<< statistics[i][j] << endl;
                send_private_batch_shares(x, sockets, nparties);
            }
        }

        cout << "Ready to send private left sums" << endl;
        vector<int> left_nums;
        for (int i = 0; i < global_split_num; i++) {
            left_nums.push_back(i+1);
            //cout << "left_nums[" << i << "] = " << left_nums[i] << endl;
        }

        for (int i = 0; i < global_split_num; i++) {
            vector<gfp> input_values_gfp(1);
            input_values_gfp[0].assign(left_nums[i]);
            send_private_inputs(input_values_gfp, sockets, nparties);
        }

        cout << "Ready to send private right sums" << endl;
        vector<int> right_nums;
        for (int i = 0; i < global_split_num; i++) {
            right_nums.push_back((i+1) * 2);
        }

        for (int i = 0; i < global_split_num; i++) {
            vector<gfp> input_values_gfp(1);
            input_values_gfp[0].assign(right_nums[i]);
            send_private_inputs(input_values_gfp, sockets, nparties);
        }

        cout << "Ready to receive result from the SPDZ engines" << endl;

        int best_split_index;
        //int best_split_index = receive_index(sockets);

        vector<float> impurities = receive_result(sockets, nparties, 3, best_split_index);
        cout << "Received: best_split_index = " << best_split_index << endl;
        cout << "Received: left_impurity = " << impurities[0] << endl;
        cout << "Received: right_impurity = " << impurities[1] << endl;

        for (unsigned int i = 0; i < sockets.size(); i++) {
            close_client_socket(sockets[i]);
        }
    }

    return 0;
}
