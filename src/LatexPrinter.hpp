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
	constexpr const char* DOCUMENT_TEMPLATE = "document.tex";

	// Outputs
	constexpr const char* TABLE_OUTPUT_DIMENSION_PREFIX = "dimension_";
	constexpr const char* LINES_TABLE_OUTPUT = "_lines_table.tex";
	constexpr const char* HYPERPLANES_TABLE_OUTPUT = "_hyperplanes_table.tex";
	constexpr const char* HYPERPLANE_REPRESENTATION_OUTPUT_PREFIX = "hyperplane_representation_";
	constexpr const char* HYPERPLANE_REPRESENTATION_OUTPUT_POSTFIX = ".tex";
	constexpr const char* DOCUMENT_OUTPUT = "tables.tex";

	// Information
	constexpr const char* TABLES_DOCUMENT_TITLE = "HyperplaneFinder tables";
}

class LatexPrinter{

public:

	LatexPrinter() noexcept
	  : m_environment(TEMPLATE_FOLDER, OUTPUT_FOLDER)
	  , m_generated_tables()
	  , m_current_hyperplane_representation_number(0){

		std::error_code ignored;
		fs::create_directories(OUTPUT_FOLDER, ignored);
		fs::create_directories(TABLES_OUTPUT_FOLDER, ignored);
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

	void generateHyperplanesTable(unsigned int geometry_dimension,
	                              const std::vector<segre::HyperplaneTableEntry>& geometry_hyp_table,
	                              size_t sub_geometries_number){
		json data;
		data["subGeometriesNumber"] = sub_geometries_number;

		std::vector<json> hyps_info;
		for(const segre::HyperplaneTableEntry& entry : geometry_hyp_table){
			json hyp_info;
			hyp_info["points"] = entry.nbrPoints;
			hyp_info["lines"] = entry.nbrLines;

			std::vector<size_t> sub;
			sub.reserve(sub_geometries_number+1);
			for(size_t i = 0; i < sub_geometries_number; ++i){
				const std::map<long long int, std::size_t>::const_iterator it = entry.subgeometry.find(static_cast<long long int>(i));
				sub.push_back(it == entry.subgeometry.end() ? 0 : it->second);
			}
			if(sub_geometries_number > 0){
				const std::map<long long int, std::size_t>::const_iterator it = entry.subgeometry.find(-1);
				sub.push_back(it == entry.subgeometry.end() ? 0 : it->second);
			}
			hyp_info["subgeometries"] = std::move(sub);

			hyp_info["cardinal"] = entry.count;

			hyps_info.push_back(std::move(hyp_info));
		}
		data["hyperplanes"] = std::move(hyps_info);

		inja::Template document = m_environment.parse_template(HYPERPLANES_TABLE_TEMPLATE);
		std::string file_path = std::string(TABLES_OUTPUT_FOLDER) + TABLE_OUTPUT_DIMENSION_PREFIX + std::to_string(geometry_dimension) + HYPERPLANES_TABLE_OUTPUT;
		m_environment.write(document, data, file_path);

		m_generated_tables.emplace_back("Dimension " + std::to_string(geometry_dimension) + " hyperplanes", file_path);
	}

	void generateTablesDocument(){
		json data;
		data["title"] = TABLES_DOCUMENT_TITLE;

		time_t now = time(nullptr);
		struct tm tstruct = *localtime(&now);
		char buf[11];
		strftime(buf, sizeof(buf), "%Y/%m/%d", &tstruct);
		data["date"] = buf;

		std::vector<json> tables;
		for(const Table& generated_table : m_generated_tables){
			json table;
			table["tableName"] = generated_table.table_name;
			table["filePath"] = generated_table.file_path;
			tables.push_back(std::move(table));
		}
		data["tables"] = std::move(tables);

		inja::Template document = m_environment.parse_template(DOCUMENT_TEMPLATE);
		m_environment.write(document, data, DOCUMENT_OUTPUT);
	}

	template<size_t Dimension, size_t NbrPointsPerLine>
	std::string generateHyperplaneRepresentation(const std::bitset<math::pow(NbrPointsPerLine,Dimension)>& hyperplane){
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

private:

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
};

#endif //HYPERPLANEFINDER_LATEXPRINTER_HPP
