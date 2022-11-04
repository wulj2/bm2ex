#include "bm2ex.h"

uint32_t opt_t::gen_cigar(uint32_t length, char op_letter){
    return (length << BAM_CIGAR_SHIFT) | (opc2opi[(int)op_letter]);
}

uint8_t* opt_t::seq2ints(const char* seq, const int len){
    uint8_t* num = (uint8_t*)malloc(len*sizeof(uint8_t));
    for(int i = 0; i < len; ++i) num[i] = nt2int[(int)seq[i]];
    return num;
}

uint32_t* opt_t::add_cigar(uint32_t* new_cigar, int32_t* p, int32_t* s, uint32_t length, char op){
    if((*p) >= (*s)){
        ++(*s);
        kroundup32(*s);
        new_cigar = (uint32_t*)realloc(new_cigar, (*s)*sizeof(uint32_t));
    }
    new_cigar[(*p)++] = gen_cigar(length, op);
    return new_cigar;
}

void opt_t::m2ex(){
    // read ref index
    faidx_t* fai = fai_load(iref.c_str());
    if(!fai){
        fprintf(stderr, "input reference must be indexed\n");
        exit(EXIT_FAILURE);
    }
    // read bam
    samFile* ifp = sam_open(ibam.c_str(), "r");
    bam_hdr_t* hdr = sam_hdr_read(ifp);
    samFile* ofp = sam_open(obam.c_str(), "wb");
    assert(sam_hdr_write(ofp, hdr) >= 0);
    uint8_t** refi = (uint8_t**)calloc(hdr->n_targets, sizeof(uint8_t*));
    int len = 0;
    char* refc = NULL;
    bam1_t* b = bam_init1();
    int32_t ltid = -1;
    while(sam_read1(ifp, hdr, b) >= 0){
        if(b->core.flag & BAM_FUNMAP){
            assert(sam_write1(ofp, hdr, b) >= 0);
        }else{
            if((b->core.tid != ltid) && refi[b->core.tid] == NULL){ // new chr
                refc = fai_fetch(fai, hdr->target_name[b->core.tid], &len);
                if(refc == NULL){
                    fprintf(stderr, "ref [%s] is not in reference [%s]\n", hdr->target_name[b->core.tid], iref.c_str());
                    exit(EXIT_FAILURE);
                }else{
                    refi[b->core.tid] = seq2ints(refc, len);
                    free(refc);
                }
                if(bs && ltid >= 0){
                    free(refi[ltid]);
                    refi[ltid] = NULL;
                }
            }
            ltid = b->core.tid;
            do1(b, refi[b->core.tid]);
            assert(sam_write1(ofp, hdr, b) >= 0);
        }
    }
    sam_close(ofp);
    sam_close(ifp);
    bam_destroy1(b);
    for(int i = 0; i < hdr->n_targets; ++i){
        if(refi[i]){
            free(refi[i]);
            refi[i] = NULL;
        }
    }
    free(refi); refi = NULL;
    sam_hdr_destroy(hdr);
    // build index if needed
    if(bs) assert(sam_index_build(obam.c_str(), 14) == 0);
}

void opt_t::do1(bam1_t* b, uint8_t* r){
    int oi, ol, oplen_x = 0, oplen_m = 0;
    int32_t s = b->core.n_cigar+2;
    int32_t p = 0;
    uint32_t* nc = (uint32_t*)malloc(s*sizeof(uint32_t));
    uint32_t* oc = bam_get_cigar(b);
    int qpos = 0;
    int rpos = b->core.pos;
    int opt = 0;
    for(uint32_t i = 0; i < b->core.n_cigar; ++i){
        oi = bam_cigar_op(oc[i]);
        ol = bam_cigar_oplen(oc[i]);
        if(oi == BAM_CMATCH){
            for(int j = 0; j < ol; ++j){
                if(seq_nt16_int[bam_seqi(bam_get_seq(b), qpos+j)] != r[rpos+j]){
                    if(oplen_m) nc = add_cigar(nc, &p, &s, oplen_m, '=');
                    oplen_m = 0;
                    ++oplen_x;
                }else{
                    if(oplen_x) nc = add_cigar(nc, &p, &s, oplen_x, 'X');
                    oplen_x = 0;
                    ++oplen_m;
                }
            }
            qpos += ol;
            rpos += ol;
        }else{
            if(oplen_m) nc = add_cigar(nc, &p, &s, oplen_m, '=');
            if(oplen_x) nc = add_cigar(nc, &p, &s, oplen_x, 'X');
            oplen_m = oplen_x = 0;
            nc = add_cigar(nc, &p, &s, ol, BAM_CIGAR_STR[oi]);
            opt = bam_cigar_type(oi);
            if(opt == 1){
                qpos += ol;
            }else if(opt == 2){
                rpos += ol;
            }
        }
    }
    if(oplen_m) nc = add_cigar(nc, &p, &s, oplen_m, '=');
    if(oplen_x) nc = add_cigar(nc, &p, &s, oplen_x, 'X');
    int32_t new_datal = b->l_data + (p<<2)-(b->core.n_cigar<<2);
    uint8_t* ndata = (uint8_t*)malloc(new_datal*sizeof(uint8_t));
    memcpy(ndata, bam_get_qname(b), b->core.l_qname);
    memcpy(ndata+b->core.l_qname, (uint8_t*)nc, p<<2);
    memcpy(ndata+b->core.l_qname+(p<<2), bam_get_seq(b), b->l_data-(b->core.l_qname+(b->core.n_cigar<<2)));
    free(b->data);
    b->data = ndata;
    b->core.n_cigar = p;
    b->m_data = new_datal;
    b->l_data = new_datal;
    free(nc);
}

void bm2ex_usage(opt_t* opt, char* arg0){
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: a lite tool to convert cigar M in BAM to =X\n\n");
    fprintf(stderr, "Usage: %s [options]\n\n", arg0);
    fprintf(stderr, "Options: -i FILE input bam\n");
    fprintf(stderr, "         -r FILE input indexed ref\n");
    fprintf(stderr, "         -o FILE output bam [%s]\n", opt->obam.c_str());
    fprintf(stderr, "         -s      input BAM is sorted by coordinates if set [no]\n");
    fprintf(stderr, "\nps: you'd better sort your input BAM by coordinates and add the -s options,\n");
    fprintf(stderr, "    this will reduce reference memory overhead,\n");
    fprintf(stderr, "    and the index for out BAM will be generated automatically for you then.\n");
    fprintf(stderr, "\n");
}

int bm2ex_main(int argc, char** argv){
    opt_t opt;
    if(argc == 1) {
        bm2ex_usage(&opt, argv[0]);
        return 0;
    }
    int c = -1;
    while((c = getopt(argc, argv, "i:r:o:sh")) >= 0) {
        switch(c) {
            case 'i': opt.ibam = optarg; break;
            case 'r': opt.iref = optarg; break;
            case 'o': opt.obam = optarg; break;
            case 's': opt.bs = true; break;
            case 'h': bm2ex_usage(&opt, argv[0]); return 0; break;
            default: break;
        }
    }
    opt.m2ex();
    return 0;
}

int main(int argc, char** argv){
    return bm2ex_main(argc, argv);
}
