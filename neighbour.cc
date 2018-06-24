
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <bitset>
#include <set>

#include <stdio.h>
#include <ctype.h>
#include <assert.h>

using Dict = std::vector<std::string>;

struct Node {
    std::vector<int> neighbours; // index of words next to this one
};

using Graph = std::vector<Node>;

using HashType = std::bitset<26> ;

static inline HashType word_hash(const std::string &a)
{
    HashType hash;
    for (auto c : a) {
        hash.set(c - 'a');
    }
    return hash;
}

static inline bool to_lower_str(std::string &s)
{
    bool good = true;
    for (auto i = s.begin(); i != s.end() && good; ++i) {
        *i = tolower(*i);
        good = islower(*i);
    }
    return good;
}

static inline int distance(const std::string &a, const std::string &b)
{
    int n = 0;
    auto i = a.begin();
    auto j = b.begin();
    while (i != a.end() && j != b.end()) {
        if (*i++ != *j++)
            ++n;
    }
    if (i != a.end()) {
        n += a.end() - i;
    } else if (j != b.end()) {
        n += b.end() - j;
    }
    return n;
}

static inline bool is_neighbour(const std::string &a, const std::string &b)
{
    return distance(a, b) == 1;
}

static inline void add_if_neighbour(const Dict &dict, unsigned i,
                                    const std::unordered_set<std::size_t> &similar_words,
                                    Graph &graph)
{
    for (std::size_t j : similar_words) {
        if (is_neighbour(dict[i], dict[j])) {
            graph[i].neighbours.push_back(j);
            graph[j].neighbours.push_back(i);
        }
    }
}

void init_graph(const Dict &dict, Graph &graph)
{
    // mapping word hashes to list of words having that hash
    using Hash2Words = std::unordered_map<HashType, std::unordered_set<std::size_t>>;

    Hash2Words hash2words;
    for (unsigned i = 0; i < dict.size(); ++i) {
        const std::string &current = dict[i];
        hash2words[word_hash(current)].insert(i);

        for (auto iter = current.begin(); iter != current.end(); ++iter) {
            std::string shorter(current.begin(), iter);
            shorter.insert(shorter.end(), iter + 1, current.end());

            hash2words[word_hash(shorter)].insert(i);
        }
    }

    // build graph
    graph.resize(dict.size());
    for (unsigned i = 0; i < dict.size(); ++i) {

        // try to find all possible neighbours by removing each character
        const std::string &current = dict[i];
        add_if_neighbour(dict, i, hash2words[word_hash(current)], graph);
        for (auto iter = current.begin(); iter != current.end(); ++iter) {
            std::string shorter(current.begin(), iter);
            shorter.insert(shorter.end(), iter + 1, current.end());

            add_if_neighbour(dict, i, hash2words[word_hash(shorter)], graph);
        }
    }
}

static std::size_t find_word(const Dict &dict, const std::string &word)
{
    auto it = std::find(dict.begin(), dict.end(), word);
    if (it == dict.end()) {
        std::cerr << "unknown word: " << word << std::endl;
        return -1;
    }
    return it - dict.begin();
}

// priority queue of nodes, shortest to the start node first
struct DistanceRecord;
using Distances = std::unordered_map<std::size_t, DistanceRecord>;

struct PriorityQueueCmp {
    PriorityQueueCmp(const Distances &distances) : distances(distances) { }
    
    inline bool operator()(std::size_t idx1, std::size_t idx2) const;

    const Distances &distances;
};

using PriorityQueue = std::multiset<std::size_t, PriorityQueueCmp>;

// mapping node indexes to shortest distance from start node,
// and its location in PriorityQueue
struct DistanceRecord {
    DistanceRecord() = default;
    
    DistanceRecord(std::size_t distance, PriorityQueue::iterator queue_ptr)
        : distance(distance), queue_ptr(queue_ptr)
    { }
    
    std::size_t	distance;
    PriorityQueue::iterator queue_ptr;
};

bool PriorityQueueCmp::operator()(std::size_t idx1, std::size_t idx2) const
{
    auto it1 = distances.find(idx1);
    auto it2 = distances.find(idx2);
    assert(it1 != distances.end());
    assert(it2 != distances.end());
    
    return it1->second.distance < it2->second.distance;
}

static bool find_path(const Dict &dict, const Graph &graph,
                      const std::string &word1, const std::string &word2)
{
    std::size_t i = find_word(dict, word1);
    if (i == -1u)
        return false;
    std::size_t j = find_word(dict, word2);
    if (j == -1u)
        return false;

    Distances distances;
    PriorityQueue queue{PriorityQueueCmp(distances)};
    auto idx_ptr = distances.emplace(i, DistanceRecord(0, queue.end())).first;
    std::size_t idx = i;
    
    while (idx_ptr != distances.end() && idx != j) {
        // add neighbours of current node 'idx'
        for (std::size_t neighbour_idx : graph[idx].neighbours) {
            auto nit = distances.find(neighbour_idx);
            if (nit == distances.end()) {
                // first time we have seen it, add to priority queue
                distances[neighbour_idx].distance = idx_ptr->second.distance + 1;
                distances[neighbour_idx].queue_ptr = queue.insert(neighbour_idx);
            } else if (nit->second.queue_ptr != queue.end()
                       && idx_ptr->second.distance + 1 < nit->second.distance) {
                // reposition in queue
                assert(nit->second.queue_ptr != queue.end());
                queue.erase(nit->second.queue_ptr);
                nit->second.distance = idx_ptr->second.distance + 1;
                nit->second.queue_ptr = queue.insert(neighbour_idx);
            }
        }

        if (!queue.empty()) {
            // take the unprocessed node with the shortest path to start
            auto qit = queue.begin();
            idx = *qit;
            queue.erase(qit);
            idx_ptr = distances.find(idx);
            idx_ptr->second.queue_ptr = queue.end();
        } else {
            idx_ptr = distances.end();
        }
    }

    return idx == j;
}

int main(int argc, char *argv[])
{
    if (argc != 4) {
        std::cerr << "usage: " << argv[0] << " <file> <word> <word>" << std::endl;
        return 1;
    }
    std::string word1(argv[2]);
    std::string word2(argv[3]);
    
    std::ifstream in(argv[1]);
    if (!in) {
        perror("open");
        return 1;
    }

    // convert to all lower case
    Dict dict;
    std::string line;
    while (std::getline(in, line)) {
        if (to_lower_str(line)) {
            dict.push_back(line);
        }
    }

    Graph graph;
    init_graph(dict, graph);

    if (!find_path(dict, graph, word1, word2)) {
        std::cerr << "not found" << std::endl;
    }

    return 0;
}
