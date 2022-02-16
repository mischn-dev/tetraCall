# code to transform diploid in tetraploid allele information

@show ARGS
if ARGS == String[]
    println(" No arguments specified \n first and only argument is the entire path to the file including the filename, e.g. ~/Documents/Variants.vvcf")

else


using CSV, DataFrames

# provide the quality filtered vcf file
vcffile = ARGS[1]
#vcffile = "/media/michael-uni/SSD_Work_Env/Nadja_test_SNPCalling/testFile.vcf"

run(`vcftools --vcf $vcffile --out diploids --max-alleles 2 --min-alleles 2 --extract-FORMAT-info AD`)

#############################################################
## if the data is not too big

vcf2 = @time CSV.read("diploids.AD.FORMAT", DataFrame; delim="\t", comment="##")

println("Data is loaded successfully, start processing")
# push each genotype in a single dict entry
call = Dict()
## the file for linking to the individual genotype dicts is processed according to the needs of CuArrays - make all strings become integers
call[:General] = Matrix(vcf2[!,1:2])
# remove the id, filter and info columns
vcf2 = vcf2[!,Not(1:2)]
GC.gc() # empty memory
genotypes = names(vcf2)
# split the individual colum for each library into 3 columns giving info about the call type, REf readcount and ALT readcount
Threads.@threads  for j in genotypes
                                  p = vcf2[!,Symbol(j)]
                                  m = Int64[]
                                  N = Int64[]
                                             for i in 1:size(p,1)
                                                                  if p[i] != "."
                                                                                push!(m,parse(Int64,string(split(p[i],",")[1]))) # refernce allele
                                                                                push!(N,parse(Int64,string(split(p[i],",")[2]))) # alternative allele
                                                                  else
                                                                      push!(m,0)
                                                                      push!(N,0)
                                                                  end
                                              end
                              call[Symbol(j)] = [m,N]
             end
vcf2 = nothing
GC.gc()
    # calulate the readdepth
    println("Calculate read depth, allele frequency and asign alleles")
    Threads.@threads     for n in genotypes
                                    p = call[Symbol(n)][1] .+ call[Symbol(n)][2]
                                    # caluclate the Af of the ref allele
                                    p2 = call[Symbol(n)][1] ./ p
                                    call[Symbol(n)] = [call[Symbol(n)][1], call[Symbol(n)][2], p, p2]

                                    # estimate the alleles for each locus
                                    alleles = String[]
                                    for b in 1:length(call[Symbol(n)][4])
                                                                            if call[Symbol(n)][4][b] >= 0.9
                                                                                                                                                push!(alleles, "0/0/0/0")
                                                                            elseif call[Symbol(n)][4][b] < 0.9 && call[Symbol(n)][4][b] > 0.625
                                                                                                                                                push!(alleles, "0/0/0/1")
                                                                            elseif call[Symbol(n)][4][b] > 0.375 && call[Symbol(n)][4][b] <= 0.625
                                                                                                                                                push!(alleles, "0/0/1/1")
                                                                            elseif call[Symbol(n)][4][b] <= 0.375 && call[Symbol(n)][4][b] > 0.1
                                                                                                                                                push!(alleles, "0/1/1/1")
                                                                            elseif call[Symbol(n)][4][b] <= 0.1
                                                                                                                                                push!(alleles, "1/1/1/1")
                                                                            else
                                                                                                                                                push!(alleles, "./././.")
                                                                            end
                                           end
                                           call[Symbol(n)] = [call[Symbol(n)][1], call[Symbol(n)][2], call[Symbol(n)][3], call[Symbol(n)][4], alleles]
              end

            ## push the allele info of all genotypes in one dataframe
            println("Create Dataframes from Vectors")
            local genofile = call[:General]
            for n in genotypes; genofile = hcat(genofile, call[Symbol(n)][5]); end
            genofile = DataFrame(genofile, :auto)
            # rename columns
            pl = ["Chr", "Pos"]
            # set the name of the snp file
            append!(pl, genotypes)
            rename!(genofile, pl)

            ## same for the allele allele frequency
            local AFfile = call[:General]
            for n in genotypes; AFfile = hcat(AFfile, call[Symbol(n)][4]); end
            AFfile = DataFrame(AFfile, :auto)
            # rename columns
            pl = ["Chr", "Pos"]
            # set the name of the snp file
            append!(pl, genotypes)
            rename!(AFfile, pl)

            ## same for the read depth
            local RDfile = call[:General]
            for n in genotypes; RDfile = hcat(RDfile, call[Symbol(n)][3]); end
            RDfile = DataFrame(RDfile, :auto)
            # rename columns
            pl = ["Chr", "Pos"]
            # set the name of the snp file
            append!(pl, genotypes)
            rename!(RDfile, pl)

            ## write output
            loc = split(vcffile, ".")[1]
            CSV.write("$(loc)_Alleles.txt", genofile; delim="\t")
            CSV.write("$(loc)_RefAF.txt", AFfile; delim="\t")
            CSV.write("$(loc)_ReadDepth.txt", RDfile; delim="\t")
            
end

rm("diploids.AD.FORMAT")
