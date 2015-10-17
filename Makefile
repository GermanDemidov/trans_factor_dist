CC=g++
CFLAGS=-c -Wall
LDFLAGS=
LPFLAGS=-pthread -std=c++11
SOURCES=main.cpp parser.cpp transcription_factor.cpp cell.cpp transcription_factors_in_cell.cpp auxuilary.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=tf_distrib

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ $(LPFLAGS)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@ $(LPFLAGS)

clean:
	rm -rf *.o








bool Cell::bind_tf_to_dna(int i, Transcription_factors_in_cell& tfs, Parser_dnase_acc& dnase) {
int choose_strand_to_bind = pick_a_number(0,1);
int coordinate = 0;
if (choose_strand_to_bind == 0) {
coordinate = pick_a_number(0, (int)forward_DNA.size() - 1);
} else {
coordinate = pick_a_number(0, (int)reverse_DNA.size() - 1);
}
std::pair<int, int> coordinates = std::make_pair(coordinate, coordinate + tfs.get_size_of_TF(i) - 1);
std::vector<double> weights_of_binding;
int size_of_tf = tfs.get_tf_by_index(i).get_size();
int new_coordinate = 0;
// maximum number of stochastic jumps before protein lands
int counter_of_stochastic_jumps = 0;

while (counter_of_stochastic_jumps < max_number_of_stochastic_jumps_before_landing) {
if (choose_strand_to_bind == 0) {
weights_of_binding = weights_of_binding_of_all_tfs_forward[(tfs.get_tf_by_index(i).type_name)];
coordinate = choose_position_with_maximum_weight_inside_neighborhood(coordinate, size_of_tf, weights_of_binding);
} else {
weights_of_binding = weights_of_binding_of_all_tfs_revcompl[(tfs.get_tf_by_index(i).type_name)];
coordinate = choose_position_with_maximum_weight_inside_neighborhood(coordinate, size_of_tf, weights_of_binding);
}

if (choose_strand_to_bind == 0) {
if (binding_site_is_free(true, coordinate, tfs.get_size_of_TF(i))) {
tfs.get_tf_by_index(i).change_coordinate_in_sequence(coordinate);
tfs.bind_tf_to_dna(i, true, coordinate);
no_access_dna_forward.insert(coordinates);
no_access_dna_revcompl.insert(reverse_interval_for_another_strand(reverse_DNA.size(), coordinates));
return true;
} else {
return false;
counter_of_stochastic_jumps++;
choose_strand_to_bind = pick_a_number(0,1);
new_coordinate = return_normal_number(coordinate, sd_of_jumps);
if (choose_strand_to_bind == 0) {
if (new_coordinate < 0 || new_coordinate > (int)forward_DNA.size() - 1) {
return false;
}
} else {
if (new_coordinate < 0 || new_coordinate > (int)reverse_DNA.size() - 1) {
return false;
}
}
if (choose_strand_to_bind != 0) {
coordinate = (int)reverse_DNA.size() - new_coordinate;
} else {
coordinate = new_coordinate;
}
}
}
else {
if (binding_site_is_free(false, coordinate, tfs.get_size_of_TF(i))) {
tfs.get_tf_by_index(i).change_coordinate_in_sequence(coordinate);
tfs.bind_tf_to_dna(i, false, coordinate);
no_access_dna_revcompl.insert(coordinates);
no_access_dna_forward.insert(reverse_interval_for_another_strand(reverse_DNA.size(), coordinates));
return true;
} else {
return false;
counter_of_stochastic_jumps++;
choose_strand_to_bind = pick_a_number(0,1);
new_coordinate = return_normal_number(coordinate, sd_of_jumps);
if (choose_strand_to_bind == 0) {
if (new_coordinate < 0 || new_coordinate > (int)forward_DNA.size() - 1) {
return false;
}
} else {
if (new_coordinate < 0 || new_coordinate > (int)reverse_DNA.size() - 1) {
return false;
}
}
if (choose_strand_to_bind != 1) {
coordinate = (int)reverse_DNA.size() - new_coordinate;
} else {
coordinate = new_coordinate;
}
}
}
}
return false;
}