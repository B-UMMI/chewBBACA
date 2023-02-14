import MUIDataTable from "mui-datatables";


const DataTable = ({ tableData, tableTitle, tableOptions, tableConditionalFormatting }) => {

	const columnData = (tableData[0].columns).map((column) => {
		return ({ name: column,
				  label: column,
				  options: {
					  filter: true,
					  sort: true,
					  display: true,
					  setCellHeaderProps: (value) => {
						  return {
							  style: {fontWeight: "bold" }
						  }
					  },
				  }
				})
	});

	// add link to locus ID
	if (tableConditionalFormatting) {
		for (let [key, value] of Object.entries(tableConditionalFormatting)) {
			columnData[key].options = {...columnData[key].options, ...value}
		}
	};

	const rowData = tableData[1].rows

	const total_data_table = (
		<MUIDataTable
			title={tableTitle}
			data={rowData}
			columns={columnData}
			options={tableOptions}
		>
		</MUIDataTable>
	  );

	return (
		<div style={{ marginTop: "40px" }}>{total_data_table}</div>
	)
};

export default DataTable;
