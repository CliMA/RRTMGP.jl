
using RRTMGP
using RRTMGP.SolarZenithAngle

@static if RRTMGP.haspkg("Plots")
  using Plots
  const export_plots = true
else
  const export_plots = false
end

function plot_day_lat_insolation(n_days::I,
                                 n_lats::I,
                                 F_arr::Array{FT},
                                 scmap,
                                 stitle,
                                 file_name) where {FT, I}
  d_arr = Array{I}(round.(linspace(0,365,num = n_days)))
  l_arr = Array{I}(round.(linspace(-90,90,num = n_lats)))
  # lats,days = np.meshgrid(l_arr,d_arr)
  # lats,days = np.meshgrid(l_arr,d_arr)

  f = 18
  # plt.figure(figsize=(15,6))
  # gs = gridspec.GridSpec(1, 2, width_ratios=[5, 1])

  # plt.subplot(gs[0])
  if scmap == "jet" || scmap == "YlOrRd"
    vmin, vmax = 0, ceil(max(F_arr...)/100)*100
  else
    mm = ceil(max(abs.(F_arr)...)/10)*10
    vmin, vmax = -mm, mm
  end

  # contourf(d_arr,l_arr,F_arr', cmap=scmap, title=stitle, vmin=vmin,vmax=vmax)
  p1 = contourf(d_arr,l_arr,F_arr', title=stitle)
  plot!(p1, size=(400,400))
  p2 = contourf(d_arr,l_arr,F_arr', title=stitle)
  plot(p1,p2, layout=(1,2), xlabel="Days since Jan 1", ylabel="Latitude")
  plot!(p1, size=(400,400))
  # contourf(d_arr,l_arr,F_arr', title=stitle, fontsize=f)
  # contourf(d_arr,l_arr,F_arr,10,cmap=scmap,vmin=vmin,vmax=vmax)
  savefig(file_name)
  # contourf!(days,lats,F_arr,10,cmap=scmap,vmin=vmin,vmax=vmax)
  # plt.contourf(days,lats,F_arr,10,cmap=scmap,vmin=vmin,vmax=vmax)

  # plt.xlabel('Days since Jan 1', fontsize=f)
  # plt.xticks(fontsize=f)
  # plt.ylabel('Latitude', fontsize=f)
  # plt.yticks(fontsize=f)
  # cbar = plt.colorbar()
  # cbar.set_label('ToA Insolation [W/m$^2$]', fontsize=f)
  # cbar.ax.tick_params(labelsize=f)
  # plt.clim([vmin,vmax])

  # plt.subplot(gs[1])
  # Fbar = sum(F_arr, dims=1)./size(F_arr,1)
  # plt.plot(Fbar,l_arr,'k-')

  # plt.xlabel('Average ToA Insolation [W/m$^2$]', fontsize=f)
  # min_insol = floor(min(mean(F_arr,axis=0))/10)*10
  # max_insol = ceil(max(mean(F_arr,axis=0))/10)*10
  # ticks = np.linspace(min_insol, max_insol, 3)
  # if max(mean(F_arr,axis=0)) < 0.5:
  #       ticks = np.linspace(-0.5, 0.5, 3)
  # plt.xticks(ticks,fontsize=f)
  # plt.yticks(visible=False)

  # plt.tight_layout()
  # plt.show()
end

