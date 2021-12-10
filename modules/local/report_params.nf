process REPORTPARAMS {
    label 'process_low'

    output:
    path(params_mqc.csv)    emit: ch_params_for_mqc

    script:
    run_params=params
    """
    echo "${run_params}" > run_params.txt
    run_params=\$(sed 's#, #\\n#g;s#^\\[##g;s#\\]\$##g;s#:#,#g' run_params.txt)

    cat <<-MQC_HEADER >> params_mqc.csv
    # section_name: 'Pipeline Parameters'
    # description: 'Final pipeline parameter configuration used in this run.'
    # format: 'csv'
    # plot_type: 'table'
    # pconfig:
    #    id: 'pipeline_parameters'
    #    table_title: "Pipeline Parameters"
    #
    parameter,value
    \$run_params
    MQC_HEADER
    """
}
