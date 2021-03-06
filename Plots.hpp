/**************************************************************************
 * Functions creating scripts for plotting with gnuplot
 **************************************************************************
 *
 * Copyright (C) Jannik Luboeinski 2020
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */


 /*** writePalViridis ***
  * Writes a palette file for gnuplot ("viridis" from Matplotlib) */
void writePalViridis()
{
	ofstream pal("palette.pal");
	pal << "# Viridis color map" << endl << endl;
	pal << "set palette defined (\\" << endl;
	pal << "0	0.267004   0.004874   0.329415,\\" << endl;
	pal << "1	0.268510   0.009605   0.335427,\\" << endl;
	pal << "2	0.269944   0.014625   0.341379,\\" << endl;
	pal << "3	0.271305   0.019942   0.347269,\\" << endl;
	pal << "4	0.272594   0.025563   0.353093,\\" << endl;
	pal << "5	0.273809   0.031497   0.358853,\\" << endl;
	pal << "6	0.274952   0.037752   0.364543,\\" << endl;
	pal << "7	0.276022   0.044167   0.370164,\\" << endl;
	pal << "8	0.277018   0.050344   0.375715,\\" << endl;
	pal << "9	0.277941   0.056324   0.381191,\\" << endl;
	pal << "10	0.278791   0.062145   0.386592,\\" << endl;
	pal << "11	0.279566   0.067836   0.391917,\\" << endl;
	pal << "12	0.280267   0.073417   0.397163,\\" << endl;
	pal << "13	0.280894   0.078907   0.402329,\\" << endl;
	pal << "14	0.281446   0.084320   0.407414,\\" << endl;
	pal << "15	0.281924   0.089666   0.412415,\\" << endl;
	pal << "16	0.282327   0.094955   0.417331,\\" << endl;
	pal << "17	0.282656   0.100196   0.422160,\\" << endl;
	pal << "18	0.282910   0.105393   0.426902,\\" << endl;
	pal << "19	0.283091   0.110553   0.431554,\\" << endl;
	pal << "20	0.283197   0.115680   0.436115,\\" << endl;
	pal << "21	0.283229   0.120777   0.440584,\\" << endl;
	pal << "22	0.283187   0.125848   0.444960,\\" << endl;
	pal << "23	0.283072   0.130895   0.449241,\\" << endl;
	pal << "24	0.282884   0.135920   0.453427,\\" << endl;
	pal << "25	0.282623   0.140926   0.457517,\\" << endl;
	pal << "26	0.282290   0.145912   0.461510,\\" << endl;
	pal << "27	0.281887   0.150881   0.465405,\\" << endl;
	pal << "28	0.281412   0.155834   0.469201,\\" << endl;
	pal << "29	0.280868   0.160771   0.472899,\\" << endl;
	pal << "30	0.280255   0.165693   0.476498,\\" << endl;
	pal << "31	0.279574   0.170599   0.479997,\\" << endl;
	pal << "32	0.278826   0.175490   0.483397,\\" << endl;
	pal << "33	0.278012   0.180367   0.486697,\\" << endl;
	pal << "34	0.277134   0.185228   0.489898,\\" << endl;
	pal << "35	0.276194   0.190074   0.493001,\\" << endl;
	pal << "36	0.275191   0.194905   0.496005,\\" << endl;
	pal << "37	0.274128   0.199721   0.498911,\\" << endl;
	pal << "38	0.273006   0.204520   0.501721,\\" << endl;
	pal << "39	0.271828   0.209303   0.504434,\\" << endl;
	pal << "40	0.270595   0.214069   0.507052,\\" << endl;
	pal << "41	0.269308   0.218818   0.509577,\\" << endl;
	pal << "42	0.267968   0.223549   0.512008,\\" << endl;
	pal << "43	0.266580   0.228262   0.514349,\\" << endl;
	pal << "44	0.265145   0.232956   0.516599,\\" << endl;
	pal << "45	0.263663   0.237631   0.518762,\\" << endl;
	pal << "46	0.262138   0.242286   0.520837,\\" << endl;
	pal << "47	0.260571   0.246922   0.522828,\\" << endl;
	pal << "48	0.258965   0.251537   0.524736,\\" << endl;
	pal << "49	0.257322   0.256130   0.526563,\\" << endl;
	pal << "50	0.255645   0.260703   0.528312,\\" << endl;
	pal << "51	0.253935   0.265254   0.529983,\\" << endl;
	pal << "52	0.252194   0.269783   0.531579,\\" << endl;
	pal << "53	0.250425   0.274290   0.533103,\\" << endl;
	pal << "54	0.248629   0.278775   0.534556,\\" << endl;
	pal << "55	0.246811   0.283237   0.535941,\\" << endl;
	pal << "56	0.244972   0.287675   0.537260,\\" << endl;
	pal << "57	0.243113   0.292092   0.538516,\\" << endl;
	pal << "58	0.241237   0.296485   0.539709,\\" << endl;
	pal << "59	0.239346   0.300855   0.540844,\\" << endl;
	pal << "60	0.237441   0.305202   0.541921,\\" << endl;
	pal << "61	0.235526   0.309527   0.542944,\\" << endl;
	pal << "62	0.233603   0.313828   0.543914,\\" << endl;
	pal << "63	0.231674   0.318106   0.544834,\\" << endl;
	pal << "64	0.229739   0.322361   0.545706,\\" << endl;
	pal << "65	0.227802   0.326594   0.546532,\\" << endl;
	pal << "66	0.225863   0.330805   0.547314,\\" << endl;
	pal << "67	0.223925   0.334994   0.548053,\\" << endl;
	pal << "68	0.221989   0.339161   0.548752,\\" << endl;
	pal << "69	0.220057   0.343307   0.549413,\\" << endl;
	pal << "70	0.218130   0.347432   0.550038,\\" << endl;
	pal << "71	0.216210   0.351535   0.550627,\\" << endl;
	pal << "72	0.214298   0.355619   0.551184,\\" << endl;
	pal << "73	0.212395   0.359683   0.551710,\\" << endl;
	pal << "74	0.210503   0.363727   0.552206,\\" << endl;
	pal << "75	0.208623   0.367752   0.552675,\\" << endl;
	pal << "76	0.206756   0.371758   0.553117,\\" << endl;
	pal << "77	0.204903   0.375746   0.553533,\\" << endl;
	pal << "78	0.203063   0.379716   0.553925,\\" << endl;
	pal << "79	0.201239   0.383670   0.554294,\\" << endl;
	pal << "80	0.199430   0.387607   0.554642,\\" << endl;
	pal << "81	0.197636   0.391528   0.554969,\\" << endl;
	pal << "82	0.195860   0.395433   0.555276,\\" << endl;
	pal << "83	0.194100   0.399323   0.555565,\\" << endl;
	pal << "84	0.192357   0.403199   0.555836,\\" << endl;
	pal << "85	0.190631   0.407061   0.556089,\\" << endl;
	pal << "86	0.188923   0.410910   0.556326,\\" << endl;
	pal << "87	0.187231   0.414746   0.556547,\\" << endl;
	pal << "88	0.185556   0.418570   0.556753,\\" << endl;
	pal << "89	0.183898   0.422383   0.556944,\\" << endl;
	pal << "90	0.182256   0.426184   0.557120,\\" << endl;
	pal << "91	0.180629   0.429975   0.557282,\\" << endl;
	pal << "92	0.179019   0.433756   0.557430,\\" << endl;
	pal << "93	0.177423   0.437527   0.557565,\\" << endl;
	pal << "94	0.175841   0.441290   0.557685,\\" << endl;
	pal << "95	0.174274   0.445044   0.557792,\\" << endl;
	pal << "96	0.172719   0.448791   0.557885,\\" << endl;
	pal << "97	0.171176   0.452530   0.557965,\\" << endl;
	pal << "98	0.169646   0.456262   0.558030,\\" << endl;
	pal << "99	0.168126   0.459988   0.558082,\\" << endl;
	pal << "100	0.166617   0.463708   0.558119,\\" << endl;
	pal << "101	0.165117   0.467423   0.558141,\\" << endl;
	pal << "102	0.163625   0.471133   0.558148,\\" << endl;
	pal << "103	0.162142   0.474838   0.558140,\\" << endl;
	pal << "104	0.160665   0.478540   0.558115,\\" << endl;
	pal << "105	0.159194   0.482237   0.558073,\\" << endl;
	pal << "106	0.157729   0.485932   0.558013,\\" << endl;
	pal << "107	0.156270   0.489624   0.557936,\\" << endl;
	pal << "108	0.154815   0.493313   0.557840,\\" << endl;
	pal << "109	0.153364   0.497000   0.557724,\\" << endl;
	pal << "110	0.151918   0.500685   0.557587,\\" << endl;
	pal << "111	0.150476   0.504369   0.557430,\\" << endl;
	pal << "112	0.149039   0.508051   0.557250,\\" << endl;
	pal << "113	0.147607   0.511733   0.557049,\\" << endl;
	pal << "114	0.146180   0.515413   0.556823,\\" << endl;
	pal << "115	0.144759   0.519093   0.556572,\\" << endl;
	pal << "116	0.143343   0.522773   0.556295,\\" << endl;
	pal << "117	0.141935   0.526453   0.555991,\\" << endl;
	pal << "118	0.140536   0.530132   0.555659,\\" << endl;
	pal << "119	0.139147   0.533812   0.555298,\\" << endl;
	pal << "120	0.137770   0.537492   0.554906,\\" << endl;
	pal << "121	0.136408   0.541173   0.554483,\\" << endl;
	pal << "122	0.135066   0.544853   0.554029,\\" << endl;
	pal << "123	0.133743   0.548535   0.553541,\\" << endl;
	pal << "124	0.132444   0.552216   0.553018,\\" << endl;
	pal << "125	0.131172   0.555899   0.552459,\\" << endl;
	pal << "126	0.129933   0.559582   0.551864,\\" << endl;
	pal << "127	0.128729   0.563265   0.551229,\\" << endl;
	pal << "128	0.127568   0.566949   0.550556,\\" << endl;
	pal << "129	0.126453   0.570633   0.549841,\\" << endl;
	pal << "130	0.125394   0.574318   0.549086,\\" << endl;
	pal << "131	0.124395   0.578002   0.548287,\\" << endl;
	pal << "132	0.123463   0.581687   0.547445,\\" << endl;
	pal << "133	0.122606   0.585371   0.546557,\\" << endl;
	pal << "134	0.121831   0.589055   0.545623,\\" << endl;
	pal << "135	0.121148   0.592739   0.544641,\\" << endl;
	pal << "136	0.120565   0.596422   0.543611,\\" << endl;
	pal << "137	0.120092   0.600104   0.542530,\\" << endl;
	pal << "138	0.119738   0.603785   0.541400,\\" << endl;
	pal << "139	0.119512   0.607464   0.540218,\\" << endl;
	pal << "140	0.119423   0.611141   0.538982,\\" << endl;
	pal << "141	0.119483   0.614817   0.537692,\\" << endl;
	pal << "142	0.119699   0.618490   0.536347,\\" << endl;
	pal << "143	0.120081   0.622161   0.534946,\\" << endl;
	pal << "144	0.120638   0.625828   0.533488,\\" << endl;
	pal << "145	0.121380   0.629492   0.531973,\\" << endl;
	pal << "146	0.122312   0.633153   0.530398,\\" << endl;
	pal << "147	0.123444   0.636809   0.528763,\\" << endl;
	pal << "148	0.124780   0.640461   0.527068,\\" << endl;
	pal << "149	0.126326   0.644107   0.525311,\\" << endl;
	pal << "150	0.128087   0.647749   0.523491,\\" << endl;
	pal << "151	0.130067   0.651384   0.521608,\\" << endl;
	pal << "152	0.132268   0.655014   0.519661,\\" << endl;
	pal << "153	0.134692   0.658636   0.517649,\\" << endl;
	pal << "154	0.137339   0.662252   0.515571,\\" << endl;
	pal << "155	0.140210   0.665859   0.513427,\\" << endl;
	pal << "156	0.143303   0.669459   0.511215,\\" << endl;
	pal << "157	0.146616   0.673050   0.508936,\\" << endl;
	pal << "158	0.150148   0.676631   0.506589,\\" << endl;
	pal << "159	0.153894   0.680203   0.504172,\\" << endl;
	pal << "160	0.157851   0.683765   0.501686,\\" << endl;
	pal << "161	0.162016   0.687316   0.499129,\\" << endl;
	pal << "162	0.166383   0.690856   0.496502,\\" << endl;
	pal << "163	0.170948   0.694384   0.493803,\\" << endl;
	pal << "164	0.175707   0.697900   0.491033,\\" << endl;
	pal << "165	0.180653   0.701402   0.488189,\\" << endl;
	pal << "166	0.185783   0.704891   0.485273,\\" << endl;
	pal << "167	0.191090   0.708366   0.482284,\\" << endl;
	pal << "168	0.196571   0.711827   0.479221,\\" << endl;
	pal << "169	0.202219   0.715272   0.476084,\\" << endl;
	pal << "170	0.208030   0.718701   0.472873,\\" << endl;
	pal << "171	0.214000   0.722114   0.469588,\\" << endl;
	pal << "172	0.220124   0.725509   0.466226,\\" << endl;
	pal << "173	0.226397   0.728888   0.462789,\\" << endl;
	pal << "174	0.232815   0.732247   0.459277,\\" << endl;
	pal << "175	0.239374   0.735588   0.455688,\\" << endl;
	pal << "176	0.246070   0.738910   0.452024,\\" << endl;
	pal << "177	0.252899   0.742211   0.448284,\\" << endl;
	pal << "178	0.259857   0.745492   0.444467,\\" << endl;
	pal << "179	0.266941   0.748751   0.440573,\\" << endl;
	pal << "180	0.274149   0.751988   0.436601,\\" << endl;
	pal << "181	0.281477   0.755203   0.432552,\\" << endl;
	pal << "182	0.288921   0.758394   0.428426,\\" << endl;
	pal << "183	0.296479   0.761561   0.424223,\\" << endl;
	pal << "184	0.304148   0.764704   0.419943,\\" << endl;
	pal << "185	0.311925   0.767822   0.415586,\\" << endl;
	pal << "186	0.319809   0.770914   0.411152,\\" << endl;
	pal << "187	0.327796   0.773980   0.406640,\\" << endl;
	pal << "188	0.335885   0.777018   0.402049,\\" << endl;
	pal << "189	0.344074   0.780029   0.397381,\\" << endl;
	pal << "190	0.352360   0.783011   0.392636,\\" << endl;
	pal << "191	0.360741   0.785964   0.387814,\\" << endl;
	pal << "192	0.369214   0.788888   0.382914,\\" << endl;
	pal << "193	0.377779   0.791781   0.377939,\\" << endl;
	pal << "194	0.386433   0.794644   0.372886,\\" << endl;
	pal << "195	0.395174   0.797475   0.367757,\\" << endl;
	pal << "196	0.404001   0.800275   0.362552,\\" << endl;
	pal << "197	0.412913   0.803041   0.357269,\\" << endl;
	pal << "198	0.421908   0.805774   0.351910,\\" << endl;
	pal << "199	0.430983   0.808473   0.346476,\\" << endl;
	pal << "200	0.440137   0.811138   0.340967,\\" << endl;
	pal << "201	0.449368   0.813768   0.335384,\\" << endl;
	pal << "202	0.458674   0.816363   0.329727,\\" << endl;
	pal << "203	0.468053   0.818921   0.323998,\\" << endl;
	pal << "204	0.477504   0.821444   0.318195,\\" << endl;
	pal << "205	0.487026   0.823929   0.312321,\\" << endl;
	pal << "206	0.496615   0.826376   0.306377,\\" << endl;
	pal << "207	0.506271   0.828786   0.300362,\\" << endl;
	pal << "208	0.515992   0.831158   0.294279,\\" << endl;
	pal << "209	0.525776   0.833491   0.288127,\\" << endl;
	pal << "210	0.535621   0.835785   0.281908,\\" << endl;
	pal << "211	0.545524   0.838039   0.275626,\\" << endl;
	pal << "212	0.555484   0.840254   0.269281,\\" << endl;
	pal << "213	0.565498   0.842430   0.262877,\\" << endl;
	pal << "214	0.575563   0.844566   0.256415,\\" << endl;
	pal << "215	0.585678   0.846661   0.249897,\\" << endl;
	pal << "216	0.595839   0.848717   0.243329,\\" << endl;
	pal << "217	0.606045   0.850733   0.236712,\\" << endl;
	pal << "218	0.616293   0.852709   0.230052,\\" << endl;
	pal << "219	0.626579   0.854645   0.223353,\\" << endl;
	pal << "220	0.636902   0.856542   0.216620,\\" << endl;
	pal << "221	0.647257   0.858400   0.209861,\\" << endl;
	pal << "222	0.657642   0.860219   0.203082,\\" << endl;
	pal << "223	0.668054   0.861999   0.196293,\\" << endl;
	pal << "224	0.678489   0.863742   0.189503,\\" << endl;
	pal << "225	0.688944   0.865448   0.182725,\\" << endl;
	pal << "226	0.699415   0.867117   0.175971,\\" << endl;
	pal << "227	0.709898   0.868751   0.169257,\\" << endl;
	pal << "228	0.720391   0.870350   0.162603,\\" << endl;
	pal << "229	0.730889   0.871916   0.156029,\\" << endl;
	pal << "230	0.741388   0.873449   0.149561,\\" << endl;
	pal << "231	0.751884   0.874951   0.143228,\\" << endl;
	pal << "232	0.762373   0.876424   0.137064,\\" << endl;
	pal << "233	0.772852   0.877868   0.131109,\\" << endl;
	pal << "234	0.783315   0.879285   0.125405,\\" << endl;
	pal << "235	0.793760   0.880678   0.120005,\\" << endl;
	pal << "236	0.804182   0.882046   0.114965,\\" << endl;
	pal << "237	0.814576   0.883393   0.110347,\\" << endl;
	pal << "238	0.824940   0.884720   0.106217,\\" << endl;
	pal << "239	0.835270   0.886029   0.102646,\\" << endl;
	pal << "240	0.845561   0.887322   0.099702,\\" << endl;
	pal << "241	0.855810   0.888601   0.097452,\\" << endl;
	pal << "242	0.866013   0.889868   0.095953,\\" << endl;
	pal << "243	0.876168   0.891125   0.095250,\\" << endl;
	pal << "244	0.886271   0.892374   0.095374,\\" << endl;
	pal << "245	0.896320   0.893616   0.096335,\\" << endl;
	pal << "246	0.906311   0.894855   0.098125,\\" << endl;
	pal << "247	0.916242   0.896091   0.100717,\\" << endl;
	pal << "248	0.926106   0.897330   0.104071,\\" << endl;
	pal << "249	0.935904   0.898570   0.108131,\\" << endl;
	pal << "250	0.945636   0.899815   0.112838,\\" << endl;
	pal << "251	0.955300   0.901065   0.118128,\\" << endl;
	pal << "252	0.964894   0.902323   0.123941,\\" << endl;
	pal << "253	0.974417   0.903590   0.130215,\\" << endl;
	pal << "254	0.983868   0.904867   0.136897,\\" << endl;
	pal << "255	0.993248   0.906157   0.143936)" << endl;

	pal.close();
}

/*** writePalRainbow ***
 * Writes a palette file for gnuplot with the color that were used for the first single channel *
 * color plots */
void writePalRainbow()
{
	ofstream pal("palette.pal");
	pal << "# Rainbow color map" << endl << endl;
	pal << "set palette defined (\\" << endl;
	pal << " 0    \"#0000FF\",\\" << endl;
	pal << " 1    \"#0000FF\",\\" << endl;
	pal << " 2    \"#0080FF\",\\" << endl;
	pal << " 3    \"#00FFFF\",\\" << endl;
	pal << " 4    \"#00FF40\",\\" << endl;
	pal << " 5    \"#80FF00\",\\" << endl;
	pal << " 6    \"#FFFF00\",\\" << endl;
	pal << " 7    \"#FFBF00\",\\" << endl;
	pal << " 8    \"#FF4000\",\\" << endl;
	pal << " 9    \"#FF0000\",\\" << endl;
	pal.close();
}

/*** createOverviewColorPlot ***
 * Creates a color plot of data that has to be specified, representing an overview over different parameter sets *
 * - dE: size of one irradiance step
 * - E_min: minimum irradiance
 * - E_max: maximum irradiance
 * - file: a partial file name for the gnuplot script and the resulting plot file
 * - ylabel: label of the y-axis
 * - ystart: lowest y-value
 * - yend: highest y-value
 * - ystep: one step of the y-axis values
 * - yformat: format (decimal places) of the y-axis tic labelling
 * - cblabel: label of the colorbar
 * - xcolumn: column in the data file in which the x values are located
 * - ycolumn: column in the data file in which the y values are located
 * - cbcolumn: column in the data file in which the color values are located
 * - threedim [optional]: three-dimensional plot (true/false) */
void createOverviewColorPlot(double dE, double E_min, double E_max, const char* file, const char *ylabel, double ystart, double yend, double ystep,
								     const char* yformat, const char *cblabel, int xcolumn, int ycolumn, int cbcolumn, bool threedim = false)
{
	ofstream gpl(concat(file, "_cplot.gpl"));
	gpl << "# Output configuration" << endl;
	gpl << "set term pdf enhanced font \"Sans, 20\" color solid lw 2.5 size 5,4" << endl;
	gpl << "set output '" << dateStr("_cplot_") << file << ".pdf'" << endl;
	gpl << "#set term png enhanced font \"Sans, 16\" size 1280,1024" << endl;
	gpl << "#set output '" << dateStr("_cplot_") << file << ".png'" << endl << endl;

	gpl << "# Color map plot configuration" << endl;
	if (!threedim)
	{
		gpl << "set view map	# use 2-dim. plot instead of a 3-dim. plot" << endl;
	}
	else
	{
		gpl << "set view 70,75 # just for 3-dim. plot to rotate" << endl;
		gpl << "set ticslevel 0 # just for 3-dim. plot to shift z=0 down into the xy plane" << endl;
	}
	gpl << "unset key" << endl;
	gpl << "load 'palette.pal'" << endl << endl;

	gpl << "# Axis configuration" << endl;
	gpl << "xmin=" << E_min << endl;
	gpl << "xmax=" << E_max << endl;
	gpl << "ymin=" << ystart << endl;
	gpl << "ymax=" << yend << endl;
	gpl << "set xtics " << dE << endl;
	gpl << "set ytics " << 2.*ystep << endl;
	gpl << "set mytics 2" << endl;
	gpl << "set xrange [xmin-" << 0.5*dE << ":xmax+" << 0.5*dE << "]" << endl; // extend for half the step size for fields in matrix
	gpl << "set yrange [ymin-" << 0.5*ystep << ":ymax+" << 0.5*ystep << "]" << endl;
	gpl << "set xlabel \"E / mW/mm²\"" << endl;
	gpl << "set ylabel \"" << ylabel << "\" offset -1,0" << endl;
	gpl << "set cblabel \"" << cblabel << "\" offset 1" << endl;
	gpl << "set format x '%.1f'" << endl;
	gpl << "set format y " << yformat << endl;
	gpl << "set format z '%.0f'" << endl << endl;

	gpl << "# Insert grid" << endl;
	gpl << "set for [i=0:((xmax-xmin)/" << dE << ")-1] arrow from (xmin+" << 0.5*dE << ")+i*" << dE << ",ymin-" << 0.5*ystep
		 << 	" to (xmin+" << 0.5*dE << ")+i*" << dE << ",ymax+" << 0.5*ystep << " nohead front lw 0.4 lc rgb \"#E1000000\"" << endl;
	gpl << "set for [i=0:((ymax-ymin)/" << ystep << ")-1] arrow from xmin-" << 0.5*dE << ",(ymin+" << 0.5*ystep << ")+i*" << ystep
		 <<   " to xmax+" << 0.5*dE << ",(ymin+" << 0.5*ystep << ")+i*" << ystep << " nohead front lw 0.4 lc rgb \"#E1000000\"" << endl << endl;


	gpl << "# Set margins" << endl;
	gpl << "set tmargin at screen 0.97" << endl;
	gpl << "set bmargin at screen 0.20" << endl;
	gpl << "set rmargin at screen 0.75" << endl;
	gpl << "set lmargin at screen 0.21" << endl << endl;

	gpl << "# Plotting" << endl;
	gpl << "set datafile missing \"nan\"" << endl;
	gpl << "splot '" << dateStr("_cplotdata.txt'") << " using " << xcolumn << ":" << ycolumn << ":" << cbcolumn << " notitle with image" << endl;
	gpl.close();

	//system(concat("gnuplot ", concat(file, "_cplot.gpl")).c_str());
}

/*** createNetworkColorPlot ***
 * Creates a network color plot of data that has to be specified *
 * - f: the gnuplot file to write to
 * - Nl: the number of neurons in one line
 * - irradiance: the irradiance (amplitude) used for the simulation (set this to a negative value if irrelevant)
 * - column: the column in the data file to be used for color-coded values
 * - prefix: the principal name of the plot
 * - postfix: if necessary, the type of stimulation (actually either "_light" or "_current")
 * - matrix: if true, a non-interpolated matrix plot is created, else an interpolated color plot
 * - cblabel: the label for the colorbar
 * - cbmin [optional]: the minimum for the colorbar
 * - cbmax [optional]: the maximum for the colorbar
 * - cbtics [optional]: one increment for the colorbar axis */
void createNetworkColorPlot(ofstream &f, int Nl, double irradiance, int column, const char* prefix, const char* postfix, bool matrix, const char* cblabel,
							double cbmin = -1.0, double cbmax = -1.0, double cbtics = -1.0)
{
	char gplcall[100];

	if (irradiance >= 0.0)
		sprintf(gplcall, "%s_%.2f", prefix, irradiance);
	else
		sprintf(gplcall, "%s", prefix);

	f << "# Output configuration" << endl;
	f << "set term pdf enhanced font \"Sans, 20\" color solid lw 2.5 size 5,4.5" << endl;
	f << "set output '" << dateStr("_") << gplcall << "_map" << postfix << ".pdf'" << endl << endl;

	f << "# Contour plot configuration" << endl;
	if (!matrix)
	{
		f << "set pm3d" << endl;
		f << "unset surface" << endl;
	}
	f << "set view map	# 2-dim. plot" << endl;
	f << "#set view 50,75 # just for 3-dim. plot to rotate # 3-dim. plot" << endl;
	f << "#set ticslevel 0 # just for 3-dim. plot to shift z=0 down into the xy plane # 3-dim. plot" << endl;
	f << "unset key" << endl;
	if (!matrix)
		f << "set pm3d interpolate 100,100 # interpolate the color (0: gnuplot automatically decides the steps to use)" << endl;
	f << "load 'palette.pal' # load palette" << endl;
	f.precision(5); // set number of significant digits in output
	if (strcmp(prefix, "firingrate") == 0)
		f << "set title \"Firing rates across the network, Ê = " << irradiance << " mW/mm²\"" << endl;
	else if (strcmp(prefix, "irradiance") == 0)
		f << "set title \"Light intensities across the network, Ê = " << irradiance << " mW/mm²\"" << endl;
	else if (strcmp(prefix, "inc_connection") == 0)
		f << "set title \"Incoming connections per neuron" << endl;
	f << "set size square" << endl << endl;  // quadratic shape

	f << "# Axis configuration" << endl;
	if (!matrix)
	{
		f << "set xrange [1" << ":" << Nl << "]" << endl;
		f << "set yrange [1" << ":" << Nl << "]" << endl;
	}
	else
	{
		f << "set xrange [0.5" << ":" << double(Nl)+0.5 << "]" << endl;
		f << "set yrange [0.5" << ":" << double(Nl)+0.5 << "]" << endl;
	}
	f << "set mxtics 2" << endl;
	f << "set mytics 2" << endl;
	if (cbmin >= 0.0 && cbmax >= 0.0)
	{
		if (cbmin == cbmax)
			f << "set cbrange [" << 0.95*cbmin << ":" << 1.05*cbmin << "]" << endl;
		else
			f << "set cbrange [" << cbmin << ":" << cbmax << "]" << endl;
	}
	if (cbtics > 0.0)
		f << "set cbtics " << cbtics << endl;
	f << "set xlabel \"Neuron column\"" << endl;
	f << "set ylabel \"Neuron row\" offset -1.5" << endl; // offset to shift the label to the left
	f << "set cblabel \"" << cblabel << "\"" << endl;
	f << "set format x '%1.0f'" << endl;
	f << "set format y '%1.0f'" << endl;
	f << "set format z '%.5f'" << endl << endl;

	f << "# Set margins" << endl;
	f << "set tmargin at screen 0.90" << endl;
	f << "set bmargin at screen 0.17" << endl;
	f << "set rmargin at screen 0.78" << endl;
	f << "set lmargin at screen 0.15" << endl << endl;

	f << "# Plotting" << endl;
	f << "set datafile missing \"nan\"" << endl;
	if (!matrix)
		f << "splot '" << dateStr("_") << gplcall << ".txt' using 1:2:" << column << " notitle with lines lt 1" << endl;
	else
		f << "splot '" << dateStr("_") << gplcall << ".txt' using 1:2:" << column << " with image" << endl;
	f.close();

	if (irradiance >= 0.0)
		sprintf(gplcall, "gnuplot %s_%.2f_map%s.gpl", prefix, irradiance, postfix);
	else
		sprintf(gplcall, "gnuplot %s_map%s.gpl", prefix, postfix);
	system(gplcall);
}


/*** createFROverStimPlot ***
 * Creates a simple plot of either the mean firing rate or the spatial standard deviation over the stimulus quantity *
 * - f: the gnuplot file to write to
 * - dE: size of one irradiance step
 * - E_max: maximum irradiance
 * - type: type of plot (types defined in the next lines, referring to columns in data file; ) */
#define PLOT_TYPE_GAUSS_AMP 2
#define PLOT_TYPE_GAUSS_SIGMA 4
#define PLOT_TYPE_MEAN 8
void createFROverStimPlot(ofstream &f, const double dE, const double E_max, const int type)
{
	f << "set term pdf enhanced font \"Sans, 20\" color solid lw 2.5" << endl;
	f << "set output '";
	if (type == PLOT_TYPE_MEAN)
		f << dateStr("_mean_fr.pdf'") << endl;
	else if (type == PLOT_TYPE_GAUSS_SIGMA)
		f << dateStr("_gauss_sigma_fr.pdf'") << endl;
	f << "set xlabel \"Ê / mW/mm^2\"" << endl;
#ifdef STIMULUS_COMPARISON
	f << "set x2label \"Î_{cst} / nA\"" << endl;
	f << "set key top left" << endl;
	f << "set key samplen 2" << endl;
	f << "set x2range [" << 0.1 << ":" << 1 << "]" << endl;
	f << "set xtics nomirror" << endl;
	f << "set x2tics" << endl;
#endif
	if (type == PLOT_TYPE_GAUSS_AMP)
		f << "set ylabel \"{/Symbol n}_{max} / Hz\"" << endl;
	else if (type == PLOT_TYPE_GAUSS_SIGMA)
		f << "set ylabel \"{/Symbol s}_{FR}\" / neuron index" << endl;
	f << "set tmargin at screen 0.95" << endl;
	f << "set bmargin at screen 0.23" << endl << endl;

	f << "plot [x=";
	if (type == PLOT_TYPE_MEAN)
		f << -0.5*dE;
	else if (type == PLOT_TYPE_GAUSS_SIGMA)
		f << +0.5*dE; // no value exists for E=0

	f << ":" << E_max+0.5*dE << "] '" << dateStr("_fr.txt'")
	  << " using 1:" << type << ":" << type+1 << " lc 1 notitle with yerrorbars, \\" << endl // data points with error bars
	  << "\t\t'" << dateStr("_fr.txt'")
	  << " using 1:" << type << " lc 1 " // lines
	  << "title \"exc.\" with lines," << " \\" << endl
	  << "\t\t'" << dateStr("_fr.txt'")
	  << " using 1:" << type+4 << ":" << type+5 << " lc 2 notitle with yerrorbars, \\" << endl // data points with error bars
	  << "\t\t'" << dateStr("_fr.txt'")
	  << " using 1:" << type+4 << " lc 2 " // lines
	  << "title \"inh.\" with lines";

#ifdef STIMULUS_COMPARISON
	f << "," << " \\" << endl
	  << "\t\t'" << dateStr("_fr.txt'")
	  << " using 12:" << type+11 << ":" << type+12 << " lc 3 notitle with yerrorbars, \\" << endl // data points with error bars
	  << "\t\t'" << dateStr("_fr.txt'")
	  << " using 12:" << type+11 << " lc 3 " // lines
	  << "title \"cst\" axes x2y1 with lines" << endl;
#endif

	f.close();
}


/*** createPeakPlot ***
 * Creates a simple plot of an average peak of either the instantaneous firing rate or the open probability over time *
 * - f: the gnuplot file to write to
 * - period_len: length of one stimulus period in ms
 * - E: irradiance for this peak plot in mW/mm²
 * - type: type of plot (types defined below by column in data file)
 * - amp: absolute amplitude of the peak
 * - mu: the "mean" of the peak (the time where the amplitude is located)
 * - base: the baseline (lowest value within the peak)
 * - t_hm_left: the time in ms at which lies the left "half maximum"
 * - t_hm_right: the time in ms at which lies the right "half maximum" */
#define PLOT_TYPE_FR 2
#define PLOT_TYPE_O 3
void createPeakPlot(ofstream &f, double period_len, double E, int type, double amp, double mu, double base, double t_hm_left, double t_hm_right)
{
	f << "set term pdf enhanced font \"Sans, 20\" color lw 2.5" << endl;

	// Set output file
	f << "set output '";
	if (type == PLOT_TYPE_FR)
		f << dateStr(string("_fr_peak_") + dtos(E,2) + string(".pdf'")) << endl << endl;
	else if (type == PLOT_TYPE_O)
		f << dateStr(string("_O_peak_") + dtos(E,2) + string(".pdf'")) << endl << endl;


	// Set variables (or rather initial values in case of Gauss fit)
	f << "A = " << amp-base << endl
	  << "mu = " << mu << endl
	  << "b = " << base << endl
	  << "FWHM = " << ((t_hm_right >= t_hm_left) ? t_hm_right-t_hm_left : period_len+t_hm_right-t_hm_left) << endl << endl;

	// Set labels
	f << "set xlabel \"t / ms\"" << endl;
	if (type == PLOT_TYPE_O)
		f << "set ylabel \"O(t)\" offset 1.8" << endl << endl;
	else if (type == PLOT_TYPE_FR)
	{
		f << "set ylabel \"{/Symbol n}_{inst} / Hz\" offset 1.8" << endl << endl;
	}
	// Print values obtained from fitting to labels and to file
	f << "# Create labels and lines" << endl
	  << "set label sprintf(\"Height = %.2f";
	if (type == PLOT_TYPE_FR)
		f << " Hz";
	f << "\", A+b) at screen 0.56,screen 0.88" << endl
	  << "set label sprintf(\"Base = %.2f";
	if (type == PLOT_TYPE_FR)
		f << " Hz";
	f << "\", b) at screen 0.56,screen 0.81 " << endl
	  << "set label sprintf(\"FWHM = %.1f ms\", FWHM) at screen 0.56,screen 0.74 " << endl;
	if (t_hm_right >= t_hm_left)
		f  << "set arrow from " << t_hm_left << ",b+A/2 to " << t_hm_right << ",b+A/2 nohead lc rgb 'blue' dt 2" << endl << endl;
	else
	{
		f  << "set arrow from " << 0 << ",b+A/2 to " << t_hm_right << ",b+A/2 nohead lc rgb 'blue' dt 2" << endl;
		f  << "set arrow from " << t_hm_left << ",b+A/2 to " << period_len << ",b+A/2 nohead lc rgb 'blue' dt 2" << endl << endl;
	}

	// Set margins
	f << "# Set margins" << endl ;
	f << "set tmargin at screen 0.95" << endl;
	f << "set bmargin at screen 0.23" << endl << endl;

	// Plot command
	f << "plot [x=0" << ":" << period_len << "] '" << dateStr("_peak_") << dtos(E,2) << ".txt' "
	  << "using 1:" << type << " notitle lc rgb 'red' with lines";

	f.close();

	// Call gnuplot script
	if (type == PLOT_TYPE_FR)
		system(concat("gnuplot fr_peak_", dtos(E,2) + string(".gpl")).c_str());
	else if (type == PLOT_TYPE_O)
		system(concat("gnuplot O_peak_", dtos(E,2) + string(".gpl")).c_str());
}

/*** createGaussFitPlot ***
 * Creates a plot of data and fit curve of a half Gaussian fit (firing rate over radius from center); *
 * a comparative gnuplot fit is also done *
 * - f: the gnuplot file to write to
 * - A: the amplitude of the Gaussian obtained by a fit
 * - sigma: the standard deviation of the Gaussian obtained by a fit
 * - b: the offset of the Gaussian obtained by a fit
 * - E: the currently adjusted irradiance
 * - center: the center of the quadratic neuron population
 * - postfix: posterior term for the filename */
void createGaussFitPlot(ofstream &f, const double A, const double sigma, const double b, const double E, const double center, string postfix)
{
	f << "set term pdf enhanced font \"Sans, 20\" color lw 1" << endl;
	f << "set output '" << dateStr("_gauss_sigma_fit_") << dtos(E,2) << postfix << ".pdf'" << endl << endl;

	f << "set xlabel \"r\"" << endl;
	f << "set ylabel \"{/Symbol n} / Hz\"" << endl;
	f << "set key samplen 1" << endl << endl;

	f << "set tmargin at screen 0.95" << endl;
	f << "set bmargin at screen 0.23" << endl << endl;

	f << "f(x) = " << A << " * exp( -x**2 / (2*" << sigma << "**2) ) + " << b << endl << endl;

	f << "g(x) = A * exp( -x**2 / (2*sigma**2) )" << endl;
	f << "A = " << A << endl;
	f << "sigma = " << sigma << endl;
	f << "fit g(x) '" << dateStr("_fr_exc_") << dtos(E,2) << ".txt' using (sqrt(($1-" << center << ")**2 + ($2-" << center << ")**2)):3:(0.2) yerror via A, sigma" << endl << endl;

	f << "plot '" << dateStr("_fr_exc_") << dtos(E,2) << ".txt' using (sqrt(($1-" << center << ")**2 + ($2-" << center << ")**2)):3:(0.2) lc 1 title 'data', \\" << endl
	  << "\t\tf(x) lc 4 lw 2 title 'GSL fit'"
	  << ", g(x) lc 2 lw 2 dt 2 title 'GPL fit'";

	f.close();

	system(concat("gnuplot gauss_sigma_fit_" + dtos(E,2) + postfix, ".gpl").c_str());
}


/*** createSRP ***
 * Creates a spike raster plot of given data *
 * - f: the gnuplot file to write to *
 * - rasterplot: a vector of the considered neuron numbers *
 * - t_max: the duration of the time series */
void createSRP(ofstream &f, const vector<int> rasterplot, const double t_max)
{
	f << "set term pdf enhanced font \"Sans, 20\" color lw 2.5" << endl;

	// Set output file
	f << "set output '" << dateStr("_spike_raster_plot.pdf'") << endl << endl;

	// Set labels
	f << "set xlabel \"t / ms\"" << endl;
	f << "unset ylabel" << endl;
	f << "set yrange [0:100]" << endl;
	f << "set ytics (";

	double ypos;
	double ypos_step;
	if (rasterplot.size() == 1)
	{
		ypos = 50.;
		ypos_step = 0.;
	}
	else
	{
		ypos = 10.; // start at y=10
		ypos_step = (90.-ypos)/(rasterplot.size()-1); // go in steps up to y=90
	}
	for (int i=0; i<rasterplot.size(); i++)
	{
		f << "\"#" << rasterplot[i] << "\" " << dtos(ypos+i*ypos_step, 2);
		if (i < rasterplot.size()-1)
			f << ", ";
	}
	f << ")" << endl << endl;

	// Plot spike time series
	f << "plot [x=0:" << t_max << "]";
	for (int i=0; i<rasterplot.size(); i++)
	{
		int col = i+3; // the third column is the first to contain spiking data
		f << "'" << dateStr("_spikes.txt") << "' using 1:($" << col << " > 0 ? " << dtos(ypos+i*ypos_step, 2) << "*$"
		  << col << " : 1/0):(\"'\") notitle with labels offset character 0,character -0.15 tc rgb '#999999'";

		if (i < rasterplot.size()-1)
			f << ", \\" << endl << "      ";
	}

	// Call gnuplot script
	f.close();
	system("gnuplot spike_raster_plot.gpl");
}
