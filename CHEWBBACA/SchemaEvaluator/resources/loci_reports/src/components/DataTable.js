import MUIDataTable from "mui-datatables";


const DataTable = ({ tableData, tableTitle, tableOptions, tableConditionalFormatting }) => {

	let columnData = (tableData[0].columns).map((column) => {
		return ({ name: column,
				  label: column,
				  options: {
					  filter: true,
					  sort: true,
					  sortThirdClickReset: true,
					  display: true,
				  },
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
		<div>{total_data_table}</div>
	)
};

export default DataTable;
