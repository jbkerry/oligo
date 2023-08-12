from oligo.design import Capture, OffTarget, Tiled

def run_capture(config, genome, fasta, bed, enzyme, oligo, blat, star_index):
    c = Capture(genome=genome, fa=fasta, config_path=config, blat=blat)
    c.gen_oligos(
        bed=bed,
        enzyme=enzyme,
        oligo=oligo,
    )
    determine_oligo_density(c, star_index)

def run_tiled(config, genome, fasta, chrom, region, contig, enzyme, step_size, oligo, blat, star_index):
    c = Tiled(genome=genome, fa=fasta, config_path=config, blat=blat)
    if contig:
        c.gen_oligos_contig(
            chrom=chrom,
            region=region,
            step=step_size,
            oligo=oligo
        )
    else:
        c.gen_oligos_capture(
            chrom=chrom,
            region=region,
            enzyme=enzyme,
            oligo=oligo
        )
    determine_oligo_density(c, star_index)

def run_off_target(config, genome, fasta, bed, step_size, max_dist, oligo, blat, star_index):
    c = OffTarget(genome=genome, fa=fasta, config_path=config, blat=blat)
    c.gen_oligos(
        bed=bed,
        step=step_size,
        max_dist=max_dist,
        oligo=oligo,
    )
    determine_oligo_density(c, star_index)

def determine_oligo_density(c, star_index):
    c.write_fasta()
    c.detect_repeats().align_to_genome(s_idx=star_index)
    c.extract_repeats().calculate_density().write_oligo_info()