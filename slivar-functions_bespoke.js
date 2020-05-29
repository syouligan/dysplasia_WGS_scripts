var config = {af_cutoff: 0.01, num_homs: 2, ch_af_cutoff: 0.1}

function hq(sample) {
if(sample.alts == -1) { return false; }
if(sample.DP < 10) { return false; }
if(sample.GQ < 10) { return false; }
if(sample.alts == 0) {
if(sample.AB > 0.02) { return false; }
return true;
}
if(sample.alts == 1) {
if(sample.AB < 0.2 || sample.AB > 0.8) { return false; }
return true;
}
if(sample.alts == 2) {
if(sample.AB < 0.98) { return false; }
return true;
}
return false;
}

function filter_pass(variant) {
if(variant.FILTER == 'PASS') { return true; }
return false;
}

function not_problematic(INFO) {
if(!('gnomadG_segdup' in INFO) && !('gnomadG_lcr' in INFO) && !('gnomadG_decoy' in INFO)) { return true; }
return false;
}

function impactful_variants(INFO) {
if(INFO.impactful) { return true; }
return false;
}

function hq_all(s) {
if(hq(s)) { return true; }
return false;
}

function hq_affected(s) {
if(s.affected) {
if(hq(s)) { return true; }
}
return !s.affected;
}

function hq_any_affected(s) {
if(s.affected) {
if(hq(s)) { return true; }
}
return false;
}

function affected_auto_dominant(s) {
if((variant.CHROM == 'X' || variant.CHROM == 'chrX') && ('gnomadG_nonpar' in INFO)) { return false; }
if (s.affected) {
return s.het;
}
return s.hom_ref;
}

function affected_auto_recessive(s) {
if((variant.CHROM == 'X' || variant.CHROM == 'chrX') && ('gnomadG_nonpar' in INFO)) { return false; }
if(s.affected){
return s.hom_alt;
}
return !s.hom_alt;
}

function affected_x_dominant(s) {
if(!(variant.CHROM == 'X' || variant.CHROM == 'chrX')) { return false; }
if(!('gnomadG_nonpar' in INFO)) { return false; }
if (s.affected) {
return s.het;
}
return s.hom_ref;
}

function affected_x_recessive(s) {
if(!(variant.CHROM == 'X' || variant.CHROM == 'chrX')) { return false; }
if(!('gnomadG_nonpar' in INFO)) { return false; }
if(s.sex == "female") {
if((s.affected && s.hom_alt) || (!s.affected && !s.hom_alt)) { return true; }
}
if(s.sex == "male") {
if((s.affected && s.het) || (s.affected && s.hom_alt) || (!s.affected && s.hom_ref)) { return true; }
}
return false;
}

function rare(s) {
if(s.hom_ref && s.affected) { return false; }
if(s.hom_alt && !s.affected) { return false; }
if(s.het && s.affected) {
if(INFO.gnomad_popmax_af <= config.af_cutoff && INFO.gnomad_nhomalt_controls <= config.num_homs) { return true; }
}
if(s.hom_alt && s.affected) {
if(INFO.gnomad_nhomalt_controls <= config.num_homs) { return true; }
}
return false;
}

function rare_ch(s) {
if((s.hom_ref || s.hom_alt) && s.affected) { return false; }
if(s.hom_alt && !s.affected) { return false; }
if(s.het && s.affected) {
if(INFO.gnomad_popmax_af <= config.ch_af_cutoff && INFO.gnomad_nhomalt_controls <= config.num_homs) { return true; }
}
return false;
}

function hq_rare(s) {
if(s.hom_ref && s.affected) { return false; }
if(s.hom_alt && !s.affected) { return false; }
if(s.affected) {
if(hq(s) && rare(s)) { return true; }
}
return false;
}

function hq_rare_ch(s) {
if((s.hom_ref || s.hom_alt) && s.affected) { return false; }
if(s.hom_alt && !s.affected) { return false; }
if(s.affected) {
if(hq(s) && rare_ch(s)) { return true; }
}
return false;
}

function hq_rare_impactful(s) {
if(s.hom_ref && s.affected) { return false; }
if(s.hom_alt && !s.affected) { return false; }
if(s.affected) {
if(hq_rare(s) && INFO.impactful) { return true; }
}
return false;
}

function hq_rare_ch_impactful(s) {
if((s.hom_ref || s.hom_alt) && s.affected) { return false; }
if(s.hom_alt && !s.affected) { return false; }
if(s.affected) {
if(hq_rare_ch(s) && INFO.impactful) { return true; }
}
return false;
}


function trio_rare(kid, INFO) {
if(kid.alts < 1) { return false; }
if(kid.alts == 1) {
if(INFO.gnomad_popmax_af <= config.af_cutoff && INFO.gnomad_nhomalt_controls <= config.num_homs) { return true; }
}
if(kid.alts == 2) {
if(INFO.gnomad_nhomalt_controls <= config.num_homs) { return true; }
}
return false;
}

function trio_hq_rare(kid, INFO) {
if(kid.alts < 1) { return false; }
if(hq(kid) && trio_rare(kid, INFO)) { return true; }
return false;
}

function trio_interesting(kid, dad, mom) {
if((kid.alts == 1 && mom.alts == 0 && dad.alts == 0) || (kid.alts == 2 && mom.alts <= 1 && dad.alts <= 1)) { return true; }
if((kid.alts == 0 || kid.alts == 2) && ((mom.alts == 2 && dad.alts == 0) || (mom.alts == 0 && dad.alts == 2))) { return true; }
return false;
}

function trio_hq_rare_interesting(kid, dad, mom, INFO) {
if(kid.alts < 1) { return false; }
if(hq(kid) && trio_rare(kid, INFO) && trio_interesting(kid, dad, mom)) { return true; }
return false;
}

function trio_hq_rare_impactful(kid, INFO) {
if(kid.alts < 1) { return false; }
if(INFO.impactful && hq(kid) && trio_rare(kid, INFO)) { return true; }
return false;
}

function trio_hq_rare_impactful_interesting(kid, dad, mom, INFO) {
if(kid.alts < 1) { return false; }
if(INFO.impactful && hq(kid) && trio_rare(kid, INFO) && trio_interesting(kid, dad, mom)) { return true; }
return false;
}

function trio_denovo(kid, dad, mom) {
if(kid.alts < 1) { return false; }
if(kid.alts >= 1 && mom.alts == 0 && dad.alts == 0) { return true; }
return false;
}

function trio_autosomal_dominant(kid, dad, mom) {
if(kid.alts < 1) { return false; }
if(kid.alts == 1 && mom.alts == 0 && dad.alts == 0) { return true; }
return false;
}

function trio_autosomal_recessive(kid, dad, mom) {
if(kid.alts < 1) { return false; }
if(kid.alts == 2 && mom.alts == 1 && dad.alts == 1) { return true; }
return false;
}

function trio_uniparent_disomy(kid, mom, dad) {
if(kid.alts == 2 && ((mom.alts == 2 && dad.alts == 0) || (mom.alts == 0 && dad.alts == 2))) { return true; }
return false;
}

function trio_x_recessive(kid, dad, mom, variant, INFO) {
if(kid.alts < 1) { return false; }
if(!(variant.CHROM == 'X' || variant.CHROM == 'chrX')) { return false; }
if(!('gnomadG_nonpar' in INFO)){ return false; }
if(kid.sex == "unknown"){ return false; }
if(mom.alts == 2 || dad.alts >= 1) { return false; }
if(kid.sex == "male") {
if(kid.alts >= 1) { return true; }
}
if(kid.sex == "female") { 
  if(kid.alts == 2) { return true; }
}
return false;
}

function trio_lenient_denovo(kid, dad, mom, INFO) {
if(kid.alts < 1) { return false; }
if(!(mom.alts == 0 && dad.alts == 0)) { return false; }
if(kid.alts >= 1 && INFO.impactful && hq(kid) && INFO.gnomad_popmax_af <= config.ch_af_cutoff && INFO.gnomad_nhomalt_controls <= config.num_homs) { return true; }
return false;
}

function trio_lenient_ar(kid, dad, mom, INFO) {
if(kid.alts < 1) { return false; }
if(mom.alts == dad.alts) { return false; }
if(mom.alts == 2 || dad.alts == 2) { return false; }
if(kid.alts >= 1 && INFO.impactful && hq(kid) && INFO.gnomad_popmax_af <= config.ch_af_cutoff && INFO.gnomad_nhomalt_controls <= config.num_homs) { return true; }
return false;
}
