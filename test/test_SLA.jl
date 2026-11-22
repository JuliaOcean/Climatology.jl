# SSH/SLA tests; need downloads, so they are skipped when _SKIP_DOWNLOADS. Included from runtests.jl.
@testset "SSH/SLA" begin

    if !_SKIP_DOWNLOADS
    SLA=read(SeaLevelAnomaly(name="sla_podaac"))
    f3=plot(SLA)
    @test isa(f3,Figure)

    file=joinpath(SLA.path,SLA.name*".nc")
    gr=SLA_PODAAC.get_grid(file=file)
    data=SLA_PODAAC.read_slice(file,gr)
    sub=SLA_PODAAC.subset(; read_from_file=file,save_to_file=true)
    @test isa(sub,String)

    SLA=read(SeaLevelAnomaly(name="sla_cmems"))
    file=joinpath(SLA.path,SLA.name*".nc")
    sub=SLA_CMEMS.subset(; read_from_file=file,save_to_file=true)
    @test isa(sub,String)
    end

end
