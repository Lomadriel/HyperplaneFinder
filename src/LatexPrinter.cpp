#include <LatexPrinter.hpp>


// Folders
const std::string LatexPrinter::Config::TEMPLATE_FOLDER = "./templates/";
const std::string LatexPrinter::Config::TABLES_TEMPLATE_FOLDER = TEMPLATE_FOLDER + "tables/";
const std::string LatexPrinter::Config::HYPERPLANES_REPRESENTATION_TEMPLATE_FOLDER = TEMPLATE_FOLDER + "representations/";

const std::string LatexPrinter::Config::OUTPUT_FOLDER = "./";
const std::string LatexPrinter::Config::TABLES_OUTPUT_FOLDER = OUTPUT_FOLDER + "tables/";
const std::string LatexPrinter::Config::HYPERPLANES_REPRESENTATIONS_OUTPUT_FOLDER = OUTPUT_FOLDER + "representations/";

// Templates
const std::string LatexPrinter::Config::LINES_TABLE_TEMPLATE = "lines_table.tex";
const std::string LatexPrinter::Config::HYPERPLANES_TABLE_TEMPLATE = "hyperplanes_table.tex";
const std::string LatexPrinter::Config::TABLES_DOCUMENT_TEMPLATE = "document.tex";

const std::string LatexPrinter::Config::HYPERPLANE_REPRESENTATION_TEMPLATE = "hyperplane_representation.tex";
const std::string LatexPrinter::Config::HYPERPLANE_REPRESENTATION_DOCUMENT_TEMPLATE = "hyperplanes_representation_document.tex";
const std::string LatexPrinter::Config::HYPERPLANE_REPRESENTATION_DIM4_PART_TEMPLATE = "hyperplane_representation_dim4_part.tex";
const std::string LatexPrinter::Config::HYPERPLANE_REPRESENTATION_DIM4_FULL_TEMPLATE = "hyperplane_representation_dim4_full.tex";

const std::string LatexPrinter::Config::HYPERPLANE_REPRESENTATION_PRINT = "hyperplane_representation_print.tex";

// Outputs
const std::string LatexPrinter::Config::TABLES_OUTPUT_PREFIX = "dimension_";
const std::string LatexPrinter::Config::LINES_TABLE_OUTPUT_POSTFIX = "_lines_table.tex";
const std::string LatexPrinter::Config::HYPERPLANES_TABLE_OUTPUT_POSTFIX = "_hyperplanes_table.tex";
const std::string LatexPrinter::Config::TABLES_DOCUMENT_OUTPUT = "tables.tex";

const std::string LatexPrinter::Config::HYPERPLANE_REPRESENTATION_OUTPUT_PREFIX = "hyperplane_representation_";
const std::string LatexPrinter::Config::HYPERPLANE_REPRESENTATION_OUTPUT_POSTFIX = ".tex";
const std::string LatexPrinter::Config::HYPERPLANE_REPRESENTATION_DIM4_OUTPUT = "hyperplane_representation.tex";
const std::string LatexPrinter::Config::HYPERPLANES_REPRESENTATION_DOCUMENT_PREFIX = "dimension_";
const std::string LatexPrinter::Config::HYPERPLANES_REPRESENTATION_DOCUMENT_OUTPUT = "_hyperplanes_representation.tex";

// Information
const std::string LatexPrinter::Config::TABLES_DOCUMENT_TITLE = "HyperplaneFinder tables";
const std::string LatexPrinter::Config::HYPERPLANES_REPRESENTATION_DOCUMENT_TITLE = "Hyperplanes representation";

// Others
const int LatexPrinter::Config::COUNT_FROM = 1;

LatexPrinter::Table::Table(const std::string& table_name_, const std::string& file_name_) noexcept
  : table_name(table_name_)
    , file_path(file_name_) {

}

LatexPrinter::LatexPrinter() noexcept
  : m_environment()
    , m_generated_tables()
    , m_current_hyperplane_representation_number(0)
    , m_current_hyperplane_representation_dim4_number(0) {

	std::error_code ignored;
	fs::create_directories(Config::OUTPUT_FOLDER, ignored);
	fs::create_directories(Config::TABLES_OUTPUT_FOLDER, ignored);
	fs::create_directories(Config::HYPERPLANES_REPRESENTATIONS_OUTPUT_FOLDER, ignored);

	m_environment.add_callback("count", 1, [this](inja::Parsed::Arguments args, json data) -> size_t {
		const size_t value = m_environment.get_argument<size_t>(args, 0, data);
		return Config::COUNT_FROM + value;
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

void LatexPrinter::generateLinesTable(unsigned int geometry_dimension,
                                      const std::vector<segre::VeldkampLineTableEntry>& geometry_lin_table,
                                      size_t points_type_number) {
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

	inja::Template document = m_environment.parse_template(
	  Config::TABLES_TEMPLATE_FOLDER
	  + Config::LINES_TABLE_TEMPLATE
	);
	std::string file_path = Config::TABLES_OUTPUT_FOLDER
	                        + Config::TABLES_OUTPUT_PREFIX
	                        + std::to_string(geometry_dimension)
	                        + Config::LINES_TABLE_OUTPUT_POSTFIX;
	m_environment.write(document, data, file_path);

	m_generated_tables.emplace_back("Dimension " + std::to_string(geometry_dimension) + " Veldkamp lines", file_path);
}
