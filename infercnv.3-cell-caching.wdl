version 1.0

workflow infercnv {
    input {
    	String output_directory
        String output_directory_stripped = sub(output_directory, "/+$", "")
        File raw_counts_matrix # the matrix of genes (rows) vs. cells (columns)
        File gene_order_file # data file containing the positions of each gene along each chromosome in the genome
        File annotations_file # a description of the cells, indicating the cell type classifications.
        String ref_group_names # the reference group names from the annotations file
        String cutoff = "1"
        Boolean denoise = false
        Boolean cluster_by_groups = false
        Boolean HMM = false
        String HMM_type = ""
        String analysis_mode = "samples"
        String tumor_subcluster_pval = "0.01"
        Boolean median_filter = false
        String additional_args = ""
        Int cpu = 4
        String memory = "32G"
        String docker = "mparikhbroad/infercnv_caching:latest"
        Int preemptible = 2
    }

    call run_infercnv_initial {
        input:
        	output_dir = output_directory_stripped,
            raw_counts_matrix = raw_counts_matrix,
            gene_order_file = gene_order_file,
            annotations_file = annotations_file,
            ref_group_names = ref_group_names,
            cutoff = cutoff,
            denoise = denoise,
            cluster_by_groups = cluster_by_groups,
            HMM = HMM,
            HMM_type = HMM_type,
            analysis_mode = analysis_mode,
            tumor_subcluster_pval = tumor_subcluster_pval,
            median_filter = median_filter,
            additional_args = additional_args,
            up_to_step = 3,
            cpu=cpu,
            memory=memory,
            docker=docker,
            preemptible=preemptible
    }

    call run_infercnv_intermediate as run_infercnv_intermediate1 {
        input:
            infercnv_tar_input = run_infercnv_initial.infercnv_tar,
        	output_dir = output_directory_stripped,
            raw_counts_matrix = raw_counts_matrix,
            gene_order_file = gene_order_file,
            annotations_file = annotations_file,
            ref_group_names = ref_group_names,
            cutoff = cutoff,
            denoise = denoise,
            cluster_by_groups = cluster_by_groups,
            HMM = HMM,
            HMM_type = HMM_type,
            analysis_mode = analysis_mode,
            tumor_subcluster_pval = tumor_subcluster_pval,
            median_filter = median_filter,
            additional_args = additional_args,
            up_to_step = 10,
            cpu=cpu,
            memory=memory,
            docker=docker,
            preemptible=preemptible
    }

    call run_infercnv_intermediate as run_infercnv_intermediate2 {
        input:
            infercnv_tar_input = run_infercnv_intermediate1.infercnv_tar,
        	output_dir = output_directory_stripped,
            raw_counts_matrix = raw_counts_matrix,
            gene_order_file = gene_order_file,
            annotations_file = annotations_file,
            ref_group_names = ref_group_names,
            cutoff = cutoff,
            denoise = denoise,
            cluster_by_groups = cluster_by_groups,
            HMM = HMM,
            HMM_type = HMM_type,
            analysis_mode = analysis_mode,
            tumor_subcluster_pval = tumor_subcluster_pval,
            median_filter = median_filter,
            additional_args = additional_args,
            up_to_step = 15,
            cpu=cpu,
            memory=memory,
            docker=docker,
            preemptible=preemptible
    }

    call run_infercnv_intermediate as run_infercnv_intermediate3 {
        input:
            infercnv_tar_input = run_infercnv_intermediate2.infercnv_tar,
        	output_dir = output_directory_stripped,
            raw_counts_matrix = raw_counts_matrix,
            gene_order_file = gene_order_file,
            annotations_file = annotations_file,
            ref_group_names = ref_group_names,
            cutoff = cutoff,
            denoise = denoise,
            cluster_by_groups = cluster_by_groups,
            HMM = HMM,
            HMM_type = HMM_type,
            analysis_mode = analysis_mode,
            tumor_subcluster_pval = tumor_subcluster_pval,
            median_filter = median_filter,
            additional_args = additional_args,
            up_to_step = 18,
            cpu=cpu,
            memory=memory,
            docker=docker,
            preemptible=preemptible
    }

    call run_infercnv_intermediate as run_infercnv_intermediate4 {
        input:
            infercnv_tar_input = run_infercnv_intermediate3.infercnv_tar,
        	output_dir = output_directory_stripped,
            raw_counts_matrix = raw_counts_matrix,
            gene_order_file = gene_order_file,
            annotations_file = annotations_file,
            ref_group_names = ref_group_names,
            cutoff = cutoff,
            denoise = denoise,
            cluster_by_groups = cluster_by_groups,
            HMM = HMM,
            HMM_type = HMM_type,
            analysis_mode = analysis_mode,
            tumor_subcluster_pval = tumor_subcluster_pval,
            median_filter = median_filter,
            additional_args = additional_args,
            up_to_step = 20,
            cpu=cpu,
            memory=memory,
            docker=docker,
            preemptible=preemptible
    }

    call run_infercnv_final {
        input:
            infercnv_tar = run_infercnv_intermediate4.infercnv_tar,
        	output_dir = output_directory_stripped,
            raw_counts_matrix = raw_counts_matrix,
            gene_order_file = gene_order_file,
            annotations_file = annotations_file,
            ref_group_names = ref_group_names,
            cutoff = cutoff,
            denoise = denoise,
            cluster_by_groups = cluster_by_groups,
            HMM = HMM,
            HMM_type = HMM_type,
            analysis_mode = analysis_mode,
            tumor_subcluster_pval = tumor_subcluster_pval,
            median_filter = median_filter,
            additional_args = additional_args,
            cpu=cpu,
            memory=memory,
            docker=docker,
            preemptible=preemptible
    }

    output {
        Array[File] infercnv_outputs = run_infercnv_final.infercnv_outputs
    }
}

task run_infercnv_initial {
    input {
        String output_dir
        File raw_counts_matrix
        File gene_order_file
        File annotations_file
        String ref_group_names
        String cutoff
        Boolean denoise
        Boolean cluster_by_groups
        Boolean HMM
        String HMM_type
        String analysis_mode
        String tumor_subcluster_pval
        Boolean median_filter = false
        String memory
        Int cpu
        String docker
        Int preemptible
        Int up_to_step
        String additional_args
    }

    command {
        set -e

        mkdir -p infercnv

        inferCNV.R \
        --raw_counts_matrix ~{raw_counts_matrix} \
        --annotations_file ~{annotations_file} \
        --gene_order_file ~{gene_order_file} \
        --num_threads ~{cpu} \
        --out_dir infercnv \
        --ref_group_names ~{ref_group_names} \
        --cutoff ~{cutoff} \
        --denoise ~{denoise} \
        --cluster_by_groups ~{cluster_by_groups} \
        --HMM ~{HMM} \
        --HMM_type ~{HMM_type} \
        --analysis_mode ~{analysis_mode} \
        --tumor_subcluster_pval ~{tumor_subcluster_pval} \
        --median_filter ~{median_filter} \
        --up_to_step ~{up_to_step} \
        ~{additional_args} \
        
        tar -zcf infercnv.tar.gz infercnv
    }

    output {
        File infercnv_tar = "infercnv.tar.gz"
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(raw_counts_matrix, "GB")*2 + 32) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}

task run_infercnv_intermediate {
    input {
    	File infercnv_tar_input
        String output_dir
        File raw_counts_matrix
        File gene_order_file
        File annotations_file
        String ref_group_names
        String cutoff
        Boolean denoise
        Boolean cluster_by_groups
        Boolean HMM
        String HMM_type
        String analysis_mode
        String tumor_subcluster_pval
        Boolean median_filter = false
        Boolean no_plot = true
        String memory
        Int cpu
        String docker
        Int preemptible
        Int up_to_step
        String additional_args
    }

    command {
        set -e

        mkdir -p infercnv

        tar -zxf ~{infercnv_tar_input} infercnv

        inferCNV.R \
        --raw_counts_matrix ~{raw_counts_matrix} \
        --annotations_file ~{annotations_file} \
        --gene_order_file ~{gene_order_file} \
        --num_threads ~{cpu} \
        --out_dir infercnv \
        --ref_group_names ~{ref_group_names} \
        --cutoff ~{cutoff} \
        --denoise ~{denoise} \
        --cluster_by_groups ~{cluster_by_groups} \
        --HMM ~{HMM} \
        --HMM_type ~{HMM_type} \
        --analysis_mode ~{analysis_mode} \
        --tumor_subcluster_pval ~{tumor_subcluster_pval} \
        --median_filter ~{median_filter} \
        --no_plot ~{no_plot}
        --up_to_step ~{up_to_step} \
        ~{additional_args} \
        
        tar -zcf infercnv.tar.gz infercnv
    }

    output {
        File infercnv_tar = "infercnv.tar.gz"
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(raw_counts_matrix, "GB")*2 + 32) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}

task run_infercnv_final {
    input {
        File infercnv_tar
        String output_dir
        File raw_counts_matrix
        File gene_order_file
        File annotations_file
        String ref_group_names
        String cutoff
        Boolean denoise
        Boolean cluster_by_groups
        Boolean HMM
        String HMM_type
        String analysis_mode
        String tumor_subcluster_pval
        Boolean median_filter = false
        String memory
        Int cpu
        String docker
        Int preemptible
        String additional_args
    }

    command {
        set -e

        mkdir -p infercnv

        tar -zxf ~{infercnv_tar} infercnv

        inferCNV.R \
        --raw_counts_matrix ~{raw_counts_matrix} \
        --annotations_file ~{annotations_file} \
        --gene_order_file ~{gene_order_file} \
        --num_threads ~{cpu} \
        --out_dir infercnv \
        --ref_group_names ~{ref_group_names} \
        --cutoff ~{cutoff} \
        --denoise ~{denoise} \
        --cluster_by_groups ~{cluster_by_groups} \
        --HMM ~{HMM} \
        --HMM_type ~{HMM_type} \
        --analysis_mode ~{analysis_mode} \
        --tumor_subcluster_pval ~{tumor_subcluster_pval} \
        --median_filter ~{median_filter} \
        ~{additional_args} \
        
        gsutil -m cp -r infercnv ${output_dir}
    }

    output {
        Array[File] infercnv_outputs = glob("infercnv/*")
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(raw_counts_matrix, "GB")*2 + 32) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}