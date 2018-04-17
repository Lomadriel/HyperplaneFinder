#ifndef HYPERPLANEFINDER_LATEXPRINTER_HPP
#define HYPERPLANEFINDER_LATEXPRINTER_HPP

#include <nlohmann/json.hpp>
#include <inja.hpp>

#include <PointGeometry.hpp>

#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

using json = nlohmann::json;

namespace{
	// Folders
	constexpr const char* TEMPLATE_FOLDER = "./templates/";

	constexpr const char* OUTPUT_FOLDER = "./";
	constexpr const char* TABLES_OUTPUT_FOLDER = "tables/";

	// Templates
	constexpr const char* LINES_TABLE_TEMPLATE = "lines_table.tex";
	constexpr const char* HYPERPLANES_TABLE_TEMPLATE = "hyperplanes_table.tex";
	constexpr const char* HYPERPLANE_REPRESENTATION_TEMPLATE = "hyperplane_representation.tex";
	constexpr const char* HYPERPLANE_REPRESENTATION_DOCUMENT_TEMPLATE = "hyperplanes_representation_document.tex";
	constexpr const char* HYPERPLANE_REPRESENTATION_DIM4_PART_TEMPLATE = "hyperplane_representation_dim4_part.tex";
	constexpr const char* HYPERPLANE_REPRESENTATION_DIM4_FULL_TEMPLATE = "hyperplane_representation_dim4_full.tex";
	constexpr const char* DOCUMENT_TEMPLATE = "document.tex";

	// Non-template files
	constexpr const char* HYPERPLANE_REPRESENTATION_PRINT = "hyperplane_representation_print.tex";

	// Outputs
	constexpr const char* TABLE_OUTPUT_DIMENSION_PREFIX = "dimension_";
	constexpr const char* LINES_TABLE_OUTPUT = "_lines_table.tex";
	constexpr const char* HYPERPLANES_TABLE_OUTPUT = "_hyperplanes_table.tex";
	constexpr const char* HYPERPLANE_REPRESENTATION_OUTPUT_PREFIX = "hyperplane_representation_";
	constexpr const char* HYPERPLANE_REPRESENTATION_OUTPUT_POSTFIX = ".tex";
	constexpr const char* HYPERPLANE_REPRESENTATION_DIM4_OUTPUT = "hyperplane_representation.tex";
	constexpr const char* HYPERPLANES_REPRESENTATION_DOCUMENT_PREFIX = "dimension_";
	constexpr const char* HYPERPLANES_REPRESENTATION_DOCUMENT_OUTPUT = "_hyperplanes_representation.tex";
	constexpr const char* DOCUMENT_OUTPUT = "tables.tex";

	// Information
	constexpr const char* TABLES_DOCUMENT_TITLE = "HyperplaneFinder tables";
	constexpr const char* HYPERPLANES_REPRESENTATION_DOCUMENT_TITLE = "Hyperplanes reprensentation";

	// Config
	constexpr int COUNT_FROM = 1;
}

template<size_t>
struct dependent_false: public std::false_type{ };

class LatexPrinter{

public:

	LatexPrinter() noexcept
	  : m_environment(TEMPLATE_FOLDER, OUTPUT_FOLDER)
	  , m_generated_tables()
	  , m_current_hyperplane_representation_number(0)
	  , m_current_hyperplane_representation_dim4_number(0) {

		std::error_code ignored;
		fs::create_directories(OUTPUT_FOLDER, ignored);
		fs::create_directories(TABLES_OUTPUT_FOLDER, ignored);

		m_environment.add_callback("count", 1, [this](inja::Parsed::Arguments args, json data) -> size_t {
			const size_t value = m_environment.get_argument<size_t>(args, 0, data);
			return COUNT_FROM + value;
		});
		m_environment.add_callback("double", 1, [this](inja::Parsed::Arguments args, json data) -> size_t {
			const size_t value = m_environment.get_argument<size_t>(args, 0, data);
			return 2 * value;
		});
		m_environment.add_callback("plusOne", 1, [this](inja::Parsed::Arguments args, json data) -> size_t {
			const size_t value = m_environment.get_argument<size_t>(args, 0, data);
			return value + 1;
		});
	}

	void generateLinesTable(unsigned int geometry_dimension,
	                        const std::vector<segre::VeldkampLineTableEntry>& geometry_lin_table,
	                        size_t points_type_number){
		json data;
		data["pointsTypeNumber"] = points_type_number;

		std::vector<json> lines_info;
		for(const segre::VeldkampLineTableEntry& entry : geometry_lin_table){
			json line_info;
			line_info["isProjective"] = entry.isProjective;
			line_info["core"]["points"] = entry.coreNbrPoints;
			line_info["core"]["lines"] = entry.coreNbrLines;

			std::vector<size_t> points_types;
			points_types.reserve(points_type_number);
			for(size_t i = 0; i < points_type_number; ++i){
				const std::map<long long int, std::size_t>::const_iterator it = entry.pointsType.find(static_cast<long long int>(i));
				points_types.push_back(it == entry.pointsType.end() ? 0 : it->second);
			}
			line_info["pointsType"] = std::move(points_types);

			line_info["cardinal"] = entry.count;

			lines_info.push_back(std::move(line_info));
		}
		data["lines"] = std::move(lines_info);

		inja::Template document = m_environment.parse_template(LINES_TABLE_TEMPLATE);
		std::string file_path = std::string(TABLES_OUTPUT_FOLDER) + TABLE_OUTPUT_DIMENSION_PREFIX + std::to_string(geometry_dimension) + LINES_TABLE_OUTPUT;
		m_environment.write(document, data, file_path);

		m_generated_tables.emplace_back("Dimension " + std::to_string(geometry_dimension) + " Veldkamp lines", file_path);
	}

	template<bool printPointsOrder, bool printSubgeometries>
	void generateHyperplanesTable(unsigned int geometry_dimension,
	                              const std::vector<segre::HyperplaneTableEntry>& geometry_hyp_table,
	                              size_t sub_geometries_number){
		json data;
		data["geometryDimension"] = geometry_dimension;
		data["ordersNumber"] = geometry_dimension + 1;
		data["subDimensionsNumber"] = geometry_dimension;
		data["subGeometriesNumber"] = sub_geometries_number;
		data["divNumber"] = geometry_dimension * (sub_geometries_number + 1);

		const bool print_subgeometries = printSubgeometries && (sub_geometries_number > 0);
		data["printSubgeometries"] = print_subgeometries;
		data["printPointsOrder"] = printPointsOrder;

		std::vector<json> hyperplanes_info;
		for(const segre::HyperplaneTableEntry& entry : geometry_hyp_table){
			json hyperplane_info;
			hyperplane_info["points"] = entry.nbrPoints;
			hyperplane_info["lines"] = entry.nbrLines;

			if(print_subgeometries){
				std::vector<std::vector<size_t>> subgeometries_by_dimensions;
				subgeometries_by_dimensions.reserve(geometry_dimension);
				for(size_t i = 0; i < geometry_dimension; ++i){
					std::vector<size_t> subgeometries;
					subgeometries.reserve(sub_geometries_number+1);
					for(size_t j = 0; j < sub_geometries_number; ++j){
						const std::map<long long int, std::size_t>::const_iterator it = entry.subgeometries[i].find(static_cast<long long int>(j));
						subgeometries.push_back(it == entry.subgeometries[i].end() ? 0 : it->second);
					}
					const std::map<long long int, std::size_t>::const_iterator it = entry.subgeometries[i].find(-1);
					subgeometries.push_back(it == entry.subgeometries[i].end() ? 0 : it->second);

					subgeometries_by_dimensions.push_back(std::move(subgeometries));
				}
				hyperplane_info["subgeometries"] = std::move(subgeometries_by_dimensions);
			}

			if constexpr(printPointsOrder){
				std::vector<unsigned int> points_order;
				points_order.reserve(geometry_dimension + 1);
				for(unsigned int i = 0; i <= geometry_dimension; ++i){
					const std::map<unsigned int, unsigned int>::const_iterator it = entry.pointsOfOrder.find(i);
					points_order.push_back(it == entry.pointsOfOrder.end() ? 0 : it->second);
				}
				hyperplane_info["pointsOrder"] = std::move(points_order);
			};

			hyperplane_info["cardinal"] = entry.count;

			hyperplanes_info.push_back(std::move(hyperplane_info));
		}
		data["hyperplanes"] = std::move(hyperplanes_info);

		inja::Template document = m_environment.parse_template(HYPERPLANES_TABLE_TEMPLATE);
		std::string file_path = std::string(TABLES_OUTPUT_FOLDER) + TABLE_OUTPUT_DIMENSION_PREFIX + std::to_string(geometry_dimension) + HYPERPLANES_TABLE_OUTPUT;
		m_environment.write(document, data, file_path);

		m_generated_tables.emplace_back("Dimension " + std::to_string(geometry_dimension) + " hyperplanes", file_path);
	}

	template<size_t Dimension, size_t NbrPointsPerLine>
	std::string generateHyperplaneRepresentation(const std::bitset<math::pow(NbrPointsPerLine,Dimension)>& hyperplane){

		if constexpr (Dimension < 4){
			return generateHyperplaneRepresentationDimensionLess4<Dimension, NbrPointsPerLine>(hyperplane);
		}
		else if constexpr (Dimension == 4){
			return generateHyperplaneRepresentationDimension4<NbrPointsPerLine>(hyperplane);
		}
		else{
			static_assert(dependent_false<Dimension>::value, "Dimension not supported");
			return std::string();
		}
	}

	template<size_t Dimension, size_t NbrPointsPerLine>
	std::string generateHyperplaneRepresentationsDocument(const std::vector<std::bitset<math::pow(NbrPointsPerLine,Dimension)>>& hyperplanes){

		json data;
		data["title"] = HYPERPLANES_REPRESENTATION_DOCUMENT_TITLE;

		time_t now = time(nullptr);
		struct tm tstruct = *localtime(&now);
		char buf[11];
		strftime(buf, sizeof(buf), "%Y/%m/%d", &tstruct);
		data["date"] = buf;

		data["dimention"] = Dimension;
		data["nbrPointsPerLine"] = NbrPointsPerLine;

		std::vector<json> hyperplanes_;
		hyperplanes_.reserve(hyperplanes.size());
		for(const std::bitset<math::pow(NbrPointsPerLine,Dimension)>& hyperplane : hyperplanes){
			json hyperplane_info;
			std::vector<unsigned int> in_points;
			in_points.reserve(hyperplane.count());
			for(unsigned int i = 0; i < math::pow(NbrPointsPerLine,Dimension); ++i){
				if(hyperplane[i]){
					in_points.push_back(i);
				}
			}
			hyperplane_info["inPoints"] = std::move(in_points);
			hyperplanes_.push_back(std::move(hyperplane_info));
		}
		data["hyperplanes"] = std::move(hyperplanes_);

		std::error_code ignored;
		fs::copy_file(
		  std::string(TEMPLATE_FOLDER) + HYPERPLANE_REPRESENTATION_PRINT,
		  std::string(OUTPUT_FOLDER) + HYPERPLANE_REPRESENTATION_PRINT,
		  ignored
		);

		inja::Template document = m_environment.parse_template(HYPERPLANE_REPRESENTATION_DOCUMENT_TEMPLATE);
		std::string output_file_name = std::string(OUTPUT_FOLDER)
		                               + HYPERPLANES_REPRESENTATION_DOCUMENT_PREFIX
		                               + std::to_string(Dimension)
		                               + HYPERPLANES_REPRESENTATION_DOCUMENT_OUTPUT;
		m_environment.write(document, data, output_file_name);

		return output_file_name;
	}

private:

	template<size_t Dimension, size_t NbrPointsPerLine>
	std::string generateHyperplaneRepresentationDimensionLess4(const std::bitset<math::pow(NbrPointsPerLine,Dimension)>& hyperplane){

		if constexpr (Dimension > 3){
			static_assert(dependent_false<Dimension>::value, "Dimension not supported");
		}

		std::vector<unsigned int> in_points;
		in_points.reserve(hyperplane.count());
		for(unsigned int i = 0; i < math::pow(NbrPointsPerLine,Dimension); ++i){
			if(hyperplane[i]){
				in_points.push_back(i);
			}
		}

		json data;
		data["inPoints"] = std::move(in_points);
		data["dimention"] = Dimension;
		data["nbrPointsPerLine"] = NbrPointsPerLine;

		inja::Template document = m_environment.parse_template(HYPERPLANE_REPRESENTATION_TEMPLATE);
		std::string output_file_name = std::string(OUTPUT_FOLDER)
		                               + HYPERPLANE_REPRESENTATION_OUTPUT_PREFIX
		                               + std::to_string(m_current_hyperplane_representation_number++)
		                               + HYPERPLANE_REPRESENTATION_OUTPUT_POSTFIX;
		m_environment.write(document, data, output_file_name);

		return output_file_name;
	}

	template<size_t NbrPointsPerLine>
	std::string generateHyperplaneRepresentationDimension4(const std::bitset<math::pow(NbrPointsPerLine,4)>& hyperplane){

		const std::string hyperplane_folder = std::string(OUTPUT_FOLDER)
		                                      + HYPERPLANE_REPRESENTATION_OUTPUT_PREFIX
		                                      + std::to_string(m_current_hyperplane_representation_dim4_number++)
		                                      + "/";
		std::error_code ignored;
		fs::create_directories(hyperplane_folder, ignored);

		std::vector<std::vector<json>> parts_infos;
		parts_infos.reserve(4);
		for(int current_dimension = 0; current_dimension < 4; ++current_dimension) {
			std::vector<json> sub_parts_infos;
			sub_parts_infos.reserve(4);
			for(int current_slice = 0; current_slice < 4; ++current_slice) {

				// Fix current dimension coord to the slice
				// non fixed dimension: 0-1-2-3
				// fixed dimension: slice
				std::vector<int> coords[4];
				for(int dimension = 0; dimension < 4; ++dimension) {
					if(dimension == current_dimension){
						coords[dimension].push_back(current_slice);
					}
					else{
						for(int slice = 0; slice < 4; ++slice) {
							coords[dimension].push_back(slice);
						}
					}
				}

				// Compute dimension multiplier
				// num = x * dimension_multiplier[0] + y * dimension_multiplier[1] + x * dimension_multiplier[2] + x * dimension_multiplier[3]
				// multiplier of the current dimension: 0
				// others dimension (with current dimension ignored): 1-4-16
				int dimension_multiplier[4];
				int current_multiply = 1;
				for(int dimension = 0; dimension < 4; ++dimension) {
					if(dimension == current_dimension){
						dimension_multiplier[dimension] = 0;
					}
					else{
						dimension_multiplier[dimension] = current_multiply;
						current_multiply *= 4;
					}
				}

				// Generate points
				std::vector<json> points;
				std::vector<json> in_points;
				int used_coords[4];
				for(size_t a = 0; a < coords[0].size(); ++a){
					used_coords[0] = coords[0][a];
					for(size_t b = 0; b < coords[1].size(); ++b){
						used_coords[1] = coords[1][b];
						for(size_t c = 0; c < coords[2].size(); ++c){
							used_coords[2] = coords[2][c];
							for(size_t d = 0; d < coords[3].size(); ++d){
								used_coords[3] = coords[3][d];

								json point;
								constexpr const char* dim_names[] = {"x", "y", "z"};
								int current_dim = 0;
								for(int i = 0; i < 4; ++i){
									if(dimension_multiplier[i] > 0){
										point[dim_names[current_dim++]] = used_coords[i];
									}
								}

								int coord_number = 0;
								for(size_t i = 0; i < 4; ++i){
									coord_number += used_coords[i] * dimension_multiplier[i];
								}
								point["coordNumber"] = coord_number;

								size_t number = 0;
								size_t multiplier = 1;
								for(size_t i = 0; i < 4; ++i){
									number += static_cast<size_t>(used_coords[i]) * multiplier;
									multiplier *= 4;
								}
								point["number"] = number;

								if(hyperplane[number]){
									in_points.push_back(point);
								}
								points.push_back(std::move(point));
							}
						}
					}
				}

				json part_infos;
				part_infos["direction"] = current_dimension;
				part_infos["slice"] = current_slice;
				part_infos["points"] = std::move(points);
				part_infos["inPoints"] = std::move(in_points);

				inja::Template document = m_environment.parse_template(HYPERPLANE_REPRESENTATION_DIM4_PART_TEMPLATE);
				const std::string output_file_name = std::string("Direction_")
				                                     + std::to_string(current_dimension)
				                                     + "_slice_"
				                                     + std::to_string(current_slice)
				                                     + ".tex";
				m_environment.write(document, part_infos, hyperplane_folder + output_file_name);
				part_infos["filePath"] = std::move(output_file_name);

				sub_parts_infos.push_back(std::move(part_infos));
			}
			parts_infos.push_back(std::move(sub_parts_infos));
		}

		json data;
		data["parts"] = parts_infos;

		inja::Template document = m_environment.parse_template(HYPERPLANE_REPRESENTATION_DIM4_FULL_TEMPLATE);
		const std::string output_file_name = hyperplane_folder + HYPERPLANE_REPRESENTATION_DIM4_OUTPUT;
		m_environment.write(document, data, output_file_name);

		return output_file_name;
	}

	struct Table{
		Table(const std::string& table_name_, const std::string& file_name_) noexcept
		  : table_name(table_name_)
		  , file_path(file_name_) {

		}

		std::string table_name;
		std::string file_path;
	};

	inja::Environment m_environment;
	std::vector<Table> m_generated_tables;

	unsigned int m_current_hyperplane_representation_number;
	unsigned int m_current_hyperplane_representation_dim4_number;
};

#endif //HYPERPLANEFINDER_LATEXPRINTER_HPP
