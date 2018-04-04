#ifndef HYPERPLANEFINDER_LATEXPRINTER_HPP
#define HYPERPLANEFINDER_LATEXPRINTER_HPP

#include <nlohmann/json.hpp>
#include <inja.hpp>

#include <PointGeometry.hpp>

using json = nlohmann::json;

namespace{
	// Folders
	constexpr const char* TEMPLATE_FOLDER = "./templates/";
	constexpr const char* OUTPUT_FOLDER = "./";

	// Templates
	constexpr const char* HYPERPLANES_TABLE_TEMPLATE = "hyperplanes_table.tex";
	constexpr const char* DOCUMENT_TEMPLATE = "document.tex";

	// Outputs
	constexpr const char* TABLE_OUTPUT_DIMENSION_PREFIX = "dimension_";
	constexpr const char* HYPERPLANES_TABLE_OUTPUT = "_hyperplanes_table.tex";
	constexpr const char* DOCUMENT_OUTPUT = "tables.tex";

	// Information
	constexpr const char* DOCUMENT_TITLE = "HyperplaneFinder tables";
}

class LatexPrinter{

public:

	LatexPrinter() noexcept
	  : m_environment(TEMPLATE_FOLDER, OUTPUT_FOLDER)
	  , m_generated_tables(){

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
		std::string file_name = TABLE_OUTPUT_DIMENSION_PREFIX + std::to_string(geometry_dimension) + HYPERPLANES_TABLE_OUTPUT;
		m_environment.write(document, data, file_name);

		m_generated_tables.emplace_back("Dimension " + std::to_string(geometry_dimension) + " hyperplanes", file_name);
	}

	void generateTablesDocument(){
		json data;
		data["title"] = DOCUMENT_TITLE;

		time_t now = time(nullptr);
		struct tm tstruct = *localtime(&now);
		char buf[11];
		strftime(buf, sizeof(buf), "%Y/%m/%d", &tstruct);
		data["date"] = buf;

		std::vector<json> tables;
		for(const Table& generated_table : m_generated_tables){
			json table;
			table["tableName"] = generated_table.tableName;
			table["fileName"] = generated_table.fileName;
			tables.push_back(std::move(table));
		}
		data["tables"] = std::move(tables);

		inja::Template document = m_environment.parse_template(DOCUMENT_TEMPLATE);
		m_environment.write(document, data, DOCUMENT_OUTPUT);
	}

private:

	struct Table{
		Table(const std::string& tableName_, const std::string& fileName_) noexcept
		  : tableName(tableName_)
		  , fileName(fileName_) {

		}

		std::string tableName;
		std::string fileName;
	};

	inja::Environment m_environment;
	std::vector<Table> m_generated_tables;
};

#endif //HYPERPLANEFINDER_LATEXPRINTER_HPP
