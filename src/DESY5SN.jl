module DESY5SN

using Artifacts, ArtifactUtils
using DelimitedFiles
using DataFrames
using LinearAlgebra
using CSV
using LaTeXStrings

function __init__()

    print("Loading data and covariance through artifact")
    global cov_data = vec(readdlm(joinpath(artifact"DESY5SN", "STAT+SYS.txt")))
    global sn_data = CSV.read(joinpath(artifact"DESY5SN", "DES-SN5YR_HD.csv"), DataFrame)

end

struct DESY5SN_data

    nSN::Integer                       # Number of SNe
    covariance::Matrix{Float64}        # Covariance matrix
    inv_covariance::Matrix{Float64}    # Inverse of covariance matrix
    data::Matrix{Float64}              # SN data matrix ("zHD","zHEL","MU","MUERR","MUERR_VPEC","MUERR_SYS")
    data_header::Vector{String}        # SN data matrix columns names ("zHD","zHEL","MU","MUERR","MUERR_VPEC","MUERR_SYS")

    function DESY5SN_data()

        if cov_data[1] != nrow(sn_data)
            throw(DimensionMismatch("Number of covariance entries does not match data entries."))
        else
            nSN = Int(cov_data[1])
        end

        data = Matrix(select(sn_data, Not(:CID, :IDSURVEY)))
        data_header = names(select(sn_data, Not(:CID, :IDSURVEY)))
        cov_mat = reshape(cov_data[2:end], nSN, nSN)
        inv_cov_mat = inv(cov_mat)
        
        new(nSN, cov_mat, inv_cov_mat, data, data_header)
    end
end

#=
#    Plotting distance modulus as function of redshift
function distance_mod_plotting(SN_data)
    dist_mod_plot = scatter(SN_data.data[:,1], SN_data.data[:,4], yerr=sqrt.(diag(SN_data.covariance)),
                markershape=:circle, markersize=4, color="dodgerblue", markerstrokecolor="navy", alpha=0.75,
                label="", xlabel=L"z", ylabel=L"\mu", title=L"\mathrm{Distance \; modulus}",
                titlefontsize=18, tickfontsize=15, guidefontsize=15,
                xformatter=:latex, yformatter=:latex)
    display(dist_mod_plot)
    return dist_mod_plot
end

#    Plotting covariance matrix
function cov_mat_plotting(SN_data) 
    covariance_plot = heatmap(log.(abs.(SN_data.covariance)), 
                  xlabel=L"\mathrm{component}", ylabel=L"\mathrm{component}", title=L"\mathrm{Log \; Abs \; Cov \; Matrix}",
                  titlefontsize=18, tickfontsize=15, guidefontsize=15
                  )#xformatter=:latex, yformatter=:latex, colorbar_formatter=:latex)
    display(covariance_plot)
    return covariance_plot
end
=#
end # module DESY5SN
